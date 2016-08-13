#!/usr/bin/env python

"""
Runs fastlmmc on Cluster Nodes
"""

import sys
import os
import subprocess
import csv
import struct
from math import ceil, log10
from pprint import pprint
import numpy as np
import logging
import textwrap
from time import time

tmpdir = 'work'
root = os.path.split(os.path.realpath(sys.argv[0]))[0]
fastlmmc = './fastlmmc'

# ignore Y chromosome
species_chroms = {'human':23, 'mouse':20}

# "submit command"
def run(cmd):
    c = subprocess.Popen(cmd, shell=True, executable='/bin/bash')
    c.communicate()


# keep and run after running fastlmm
def process_output(input_fn, output_fn, numeric, igv):

    prefix, ext = os.path.splitext(os.path.basename(input_fn))

    selected_fields = ['SNP', 'Pvalue', 'SNPWeight', 'OddsRatio']
    igv_fields = ['SNP', 'Chromosome', 'Position', 'Pvalue', 'SNPWeight']


    if species == 'mouse': # snp IDs are constant 11-character strings
        # format strings must match the fields specified above
        if numeric:
            # use numeric index directly
            probe_id = int(prefix)
            #fmt = '<i11sdddd'
            fmt = '<i11sddd'
        else:
            probe_id = prefix.replace('^','/')
            # probe_id will be constant for a given file; snp_id is always 11 characters
            #fmt = '<h%ss11sdddd' % len(probe_id)
            fmt = '<h%ss11sddd' % len(probe_id)

        converter = struct.Struct(fmt)

    else: # snp IDs are rsIDs (need to adjust the format strings below); generate packer for every length between 5 and 11
        converter = None
        converters = {}

        if numeric:
            # use numeric index directly
            probe_id = int(prefix)
            for i in range(4, 12):
                fmt = '<ih%ssddd' % i
                converters[i] = struct.Struct(fmt)
        else:
            probe_id = prefix.replace('^','/')
            for i in range(4, 12):
                fmt = '<h%ssh%ssddd' % (len(probe_id), i)
                converters[i] = struct.Struct(fmt)


    f = open(input_fn)

    if igv:
        igv_file = open(output_fn.replace('.bcp', '.gwas'), 'w')
        igv_file.write('SNP\tCHR\tBP\tP\tBeta\n')
        igv_lines = []

    d = csv.DictReader(f, dialect='excel-tab')
    tmp = []
    lines = []

    for row in d:
        if igv:
            igv_lines.append([row[k] for k in igv_fields])
        else:
            if numeric:
                lines.append([probe_id] + [row[k] for k in selected_fields])
            else:
                lines.append([len(probe_id), probe_id] + [row[k] for k in selected_fields])

    # gwas file
    if igv:
        igv_lines.sort(key=lambda x:(int(x[1]), int(x[2])))
        for row in igv_lines:
            # this works for human and mouse; haven't checked other species yet
            if row[1] == str(species_chroms[species]):
                row[1] == 'X'

        igv_file.write('\n'.join(['\t'.join(row) for row in igv_lines]))
        igv_file.close()

    # bcp file
    else:
        lines.sort(key=lambda s:s[2-(numeric > 0)].lower())
        for line in lines:
            if species == 'mouse': # fixed-length JAX snp IDs
                #line = line[:-4] + map(float, line[-4:])
                line = line[:-3] + list( map(float, line[-3:]) )
                # kludge for SQL Server (infinity not supported)
                for i, item in enumerate(line):
                    if item == float('INF'):
                        line[i] = 1.79e+308
                # tmp.append(converter.pack(*line)) TODO - need sql server?

            else:
                # get length of rsID, use appropriate packer
                n = len(line[-4])
                line = line[:-4] + [n, line[-4]] + list( map(float, line[-3:]) )
                for i, item in enumerate(line):
                    if item == float('INF'):
                        line[i] = 1.79e+308
                # tmp.append(converters[n].pack(*line))


        bcp = open(output_fn, 'wb')
        # bcp.write(''.join(tmp))
        bcp.close()

    f.close()


# run fastlmmc
def run_fastlmmc(dataset, output_dir, pheno_index, covFile=None, numeric=0, igv=False, species='mouse', maxthreads=1, keeptxt=False, featsel=False, exclude=False, condition=None):
    fastlmm_cmd = textwrap.dedent('''\
        %(fastlmmc)s -bfile %(bfile)s -bfileSim %(bfile)s
        -pheno %(pheno_file)s -mpheno %(pheno_index)s
        -setOutputPrecision 5 -out %(output_fn)s -v
        -extract %(include_file)s -extractSim %(exclude_file)s
        %(condition)s
        -maxThreads %(maxthreads)s %(covar)s''').replace('\n', ' ')

    if exclude:
        fastlmm_cmd = textwrap.dedent('''\
            %(fastlmmc)s -bfile %(bfile)s -bfileSim %(bfile)s
            -pheno %(pheno_file)s -mpheno %(pheno_index)s
            -setOutputPrecision 5 -out %(output_fn)s -v
            -extract %(include_file)s -extractSim %(exclude_file)s
            %(condition)s
            -maxThreads %(maxthreads)s %(covar)s -excludeByPosition 2000000
            ''').replace('\n', ' ')

    bfile = dataset
    pheno_file = '%s.pheno.txt' % dataset

    if covFile:
        covar = '-covar %s' % covFile
    else:
        covar = ''

    f = open(pheno_file)
    # skip FID and IID columns
    headers = f.readline().strip().split('\t')[2:]
    f.close()

    # numeric should be 0, 1, or 2; enforced by fastlmm_pipeline.py
    # 0 -> use phenotype name
    # 1 -> use automatic numbering
    # 2 -> assume phenotype name is numeric and encode it as such
    if numeric == 1:
        n_digits = int(ceil(log10(len(headers))))
        fmt = '%%0%sd' % n_digits
        phenotype_name = fmt % pheno_index
    elif numeric == 2:
        numeric_headers = map(int, headers)
        n_digits = int(ceil(log10(max(numeric_headers))))
        fmt = '%%0%sd' % n_digits
        phenotype_name = fmt % numeric_headers[pheno_index - 1]
    else:
        phenotype_name = headers[pheno_index - 1].replace('/', '^')

    v = globals()

    chroms = map(str, range(1, species_chroms[species] + 1))
    if condition:
        condition = '-SnpId1 %s' % condition[0]
    else:
        condition = ''

    # temporary kludge because -excludeByPosition option is slow (at least for v2.05 and v2.06)

    v.update(locals())
    if not os.path.exists('%(output_dir)s/%(phenotype_name)s.bcp' % v) and not os.path.exists('%(output_dir)s/%(phenotype_name)s.gwas' % v):
        for i, chrom in enumerate(chroms):
            include_file = os.path.join(bfile, 'filtered.snp_ids.chr.%s.included.txt' % chrom)
            if featsel:
                exclude_file = os.path.join(tmpdir, '%s.selected_snps.chr.%s.excluded.txt' % (pheno_index, chrom))
            else:
                exclude_file = os.path.join(bfile, 'filtered.snp_ids.chr.%s.excluded.txt' % chrom)
            output_fn = os.path.join(tmpdir, '%s.%02d.fastlmm.txt' % (phenotype_name, i+1))
            v.update(locals())
            cmd = fastlmm_cmd % v
            print('\n%s' % cmd)
            sys.stdout.flush()
            subprocess.check_call(cmd, shell=True)

        # merge txt files: keep header line from first file, skip it for all others
        merge_txt_cmd = 'cat <(ls %(tmpdir)s/%(phenotype_name)s.*.txt | head -1 | xargs head -n1) <(tail -n+2 -q %(tmpdir)s/%(phenotype_name)s.*.txt) > %(output_dir)s/%(phenotype_name)s.txt' % v
        print(merge_txt_cmd)
        run(merge_txt_cmd)

        # convert merged txt file to binary format for loading into SQL Server / generate .GWAS file for loading into IGV
        process_output('%(output_dir)s/%(phenotype_name)s.txt' % v,
                       '%(output_dir)s/%(phenotype_name)s.bcp' % v,
                       numeric,
                       igv)

        if not keeptxt:
            # remove merged text file to save space
            os.remove('%(output_dir)s/%(phenotype_name)s.txt' % v)



if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    #parser.set_usage('''%prog [options] dataset''')
    parser.add_argument('-s', '--species', dest='species', help='mouse or human',
                        default=None, action='store')
    parser.add_argument('-c', '--covariate_file', dest='covFile', help='use covariate file',
                        default=None, action='store')
    parser.add_argument('-n', '--numeric', dest='numeric', help='use numeric phenotype IDs',
                        nargs='?', default=0, const=1, choices=[0, 1, 2], type=int, action='store')
    parser.add_argument('-i', '--igv-output', dest='igv', help='create IGV tracks',
                        default=False, action='store_true')
    parser.add_argument('-f', '--feature-selection', dest='featsel', help='perform feature selection',
                        default=False, action='store_true')
    parser.add_argument('-e', '--excludeByPosition', dest='exclude', help='exclude SNPs within 2Mb of tested SNP from kinship matrix construction',
                        default=False, action='store_true')
    parser.add_argument('--maxthreads', dest='maxthreads', help='max # of threads to use',
                        default=1, choices=range(1,17), type=int, action='store')
    parser.add_argument('--keeptxt', dest='keeptxt', help='keep merged text file for each phenotype',
                        default=False, action='store_true')
    parser.add_argument('--condition', dest='condition', help='condition on a SNP',
                        default=None, action='store', nargs=1)
    parser.add_argument('--debug', dest='debug', help='log debugging output',
                        default=False, action='store_true')
    parser.add_argument('dataset', help='dataset to run', action='store')
    parser.add_argument('pheno_index', help= 'phenotype index', action='store')

    args = parser.parse_args()

    species = args.species
    covFile = args.covFile
    numeric = args.numeric
    igv = args.igv
    featsel = args.featsel
    exclude = args.exclude
    maxthreads = args.maxthreads
    keeptxt = args.keeptxt
    condition = args.condition
    debug = args.debug
    dataset = args.dataset
    output_dir = root
    pheno_index = int( args.pheno_index ) + 1

    if debug:
        print('args:')
        pprint(args)

    run_fastlmmc(dataset, output_dir, pheno_index, covFile, numeric, igv, species, maxthreads, keeptxt=keeptxt, featsel=featsel, exclude=exclude, condition=condition)



# for when they fix the slow -excludeByPosition performance
#output_fn = os.path.join(output_dir, '%s.fastlmm.txt' % phenotype_name.replace('/', '^'))

#cmd = fastlmm_cmd % (globals())
#c = subprocess.Popen(cmd, shell=True)
#c.communicate()

#process_output(output_fn)
