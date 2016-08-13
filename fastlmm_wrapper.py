#!/usr/bin/env python

import sys
import os
import subprocess
import csv
import struct
from math import ceil, log10
from pprint import pprint
import logging
import textwrap
from time import time
import numpy as np
from fastlmm.association import single_snp
import pysnptools.util
from pysnptools.util.pheno import loadOnePhen
from pysnptools.snpreader import Bed

root = os.path.split(os.path.realpath(sys.argv[0]))[0]
fastlmmc = './fastlmmc'

# ignore Y chromosome
species_chroms = {'human':23, 'mouse':20}


# run fastlmmc
def run_fastlmmc(dataset, output_dir, pheno_index, covFile=None, numeric=0, igv=False, species='mouse', maxthreads=1, featsel=False, exclude=False, condition=None):
    
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
    snp_reader = Bed(bfile)

    pheno_file = loadOnePhen('%s.pheno.txt' % dataset, i_pheno = pheno_index)
    phenotype_name = pheno_file['header'][0]

    if covFile:
        covar = '-covar %s' % covFile
    else:
        covar = ''

    # numeric should be 0, 1, or 2; enforced by fastlmm_pipeline.py
    # 0 -> use phenotype name
    # 1 -> use automatic numbering
    # 2 -> assume phenotype name is numeric and encode it as such
    # if numeric == 1:
    #     n_digits = int(ceil(log10(len(headers))))
    #     fmt = '%%0%sd' % n_digits
    #     phenotype_name = fmt % pheno_index
    # elif numeric == 2:
    #     numeric_headers = map(int, headers)
    #     n_digits = int(ceil(log10(max(numeric_headers))))
    #     fmt = '%%0%sd' % n_digits
    #     phenotype_name = fmt % numeric_headers[pheno_index - 1]
    # else:
    #     phenotype_name = headers[pheno_index - 1].replace('/', '^')

    v = globals()

    chroms = map(str, range(1, species_chroms[species] + 1))
    
    # if condition:
    #     condition = '-SnpId1 %s' % condition[0]
    # else:
    #     condition = ''

    # temporary kludge because -excludeByPosition option is slow (at least for v2.05 and v2.06)

    v.update(locals())
    
    # loop through chromosomes and run
    if not os.path.exists('%(output_dir)s/%(phenotype_name)s.bcp' % v) and not os.path.exists('%(output_dir)s/%(phenotype_name)s.gwas' % v):
        for i, chrom in enumerate(chroms):
            
        # separate by chromosome for LOOCV
        test_snps = snp_reader[:, snp_reader.pos[:,0] == int(chrom)]
        mat_snps = snp_reader[:, snp_reader.pos[:,0] != int(chrom)]

        # run snp with covar
        if covFile:
            df = single_snp(test_snps = test_snps, pheno = pheno_file, K0 = mat_snps, covar = covFile)
        else:
            df = single_snp(test_snps = test_snps, pheno = pheno_file, K0 = mat_snps)

        # format outputs
        out_df = df.loc[:, ['SNP', 'Chr', 'ChrPos', 'PValue', 'SnpWeight']]
        out_df.columns = ['SNP', 'CHR', 'BP', 'P', 'Beta']
	    
        # save results into data frame
        if i == 0:
            final = out_df
        else:
            final = final.append(out_df)
                

    # output to csv
    final.to_csv('%(output_dir)s/%(phenotype_name)s.gwas' % v, sep='\t', index=False)


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
    condition = args.condition
    debug = args.debug
    dataset = args.dataset
    output_dir = root
    pheno_index = int( args.pheno_index ) 

    if debug:
        print('args:')
        pprint(args)

    run_fastlmmc(dataset, output_dir, pheno_index, covFile, numeric, igv, species, maxthreads, featsel=featsel, exclude=exclude, condition=condition)


