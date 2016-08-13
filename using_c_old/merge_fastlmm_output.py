#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import stat
import subprocess

root = os.path.split(os.path.split(os.path.realpath(sys.argv[0]))[0])[0]
prog_path = os.path.join(root, 'scripts')
dataLoc = os.path.join(root, 'data')

def process(folder, numeric=0, memory=1200, check=False, igv=False, species='mouse', dryrun=False):
    print('Processing %s' % folder)
    folder = folder.rstrip(os.path.sep)
    folder_parent, folder_base = os.path.split(folder)

    pheno_key = open(os.path.join(dataLoc, folder_base, 'pheno.key.txt'))
    phenos = dict([x.strip().split('\t')[::-1] for x in pheno_key.readlines()])
    pheno_key.close()

    if igv:
        gwas = set([x.replace('.gwas', '') for x in os.listdir(folder) if x.endswith('.gwas')])

        if numeric:
            missing = set(phenos.values()).difference(bcps)
        else:
            missing_tmp = set(phenos.keys()).difference(gwas)
            missing = set([phenos[k] for k in missing_tmp])

    else:

        bcps = set([x.replace('.bcp', '') for x in os.listdir(folder) if x.endswith('.bcp')])
        if numeric == 0:
            bcps = [x.replace('^', '/') for x in bcps]
            numeric_bcp_len = 0
        else:
            # expected row size for .bcp output file
            # must match the output format in fastlmm_wrapper.py
            # currently 4-byte integer + 11-character SNP ID + 3 float values (8 bytes each)
            #  need to modify for non-constant SNP ID lengths
            numeric_bcp_len = 4 + 11 + 24

        id_sizes = {}
        for pheno in phenos:
            try:
                id_sizes[len(pheno)] += 1
            except KeyError:
                id_sizes[len(pheno)] = 1

        row_sizes = dict([(n, 2+n+11+24) for n in id_sizes])

        # Retrieve # of SNPs that passed filtering critera
        f = open('%s/%s.bim' % (dataLoc, folder_base))
        n_snps = len(f.readlines())
        f.close()

        print('Calculating expected file sizes')
        if numeric:
            expected_sizes = dict(((bcp, numeric_bcp_len * n_snps) for bcp in bcps))
        else:
            expected_sizes = dict(((pheno, row_sizes[len(pheno)] * n_snps) for pheno in phenos))

        actual_sizes = {}
        print('Getting actual file sizes')
        for i, bcp in enumerate(bcps):
            # show a period for each percent complete
            if len(bcps) > 1000 and i % (len(bcps) / 100) == 0 and i:
                sys.stdout.write('.')
                sys.stdout.flush()
            if numeric:
                actual_sizes[bcp] = os.stat(os.path.join(folder, '%s.bcp' % bcp))[stat.ST_SIZE]
            else:
                actual_sizes[bcp] = os.stat(os.path.join(folder, '%s.bcp' % bcp.replace('/', '^')))[stat.ST_SIZE]

        print
        sys.stdout.flush()

        if len(bcps) == len(phenos) and (all((actual_sizes[bcp] == expected_sizes[bcp] for bcp in bcps)) or check is False):
            # all expected output files are present and have the right sizes
            print('Merging output files')
            cmd = 'cat %s/*.bcp > %s/%s.bcp' % (folder, folder_parent, folder_base)
            print(cmd)
            sys.stdout.flush()
            e = subprocess.call(cmd, shell=True)
        else:
            # need to re-run some tasks
            #  Most likely cause of tasks not finishing properly (assuming properly formatted input)
            #  is insufficient memory requested
            #  Therefore, double amount requested last time, unless it exceeds 4096
            if memory == 4096:
                sys.stderr.write('Previous attempt already used 4GB; some other issue may be at fault')
                sys.exit(-1)

            memory = min(memory * 2, 44000)

            if numeric:
                missing = set(phenos.values()).difference(bcps)
                wrongsize = set([bcp for bcp in bcps if actual_sizes[bcp] != expected_sizes[bcp]])
                if check is False:
                    wrongsize = set()
                sys.stderr.write('Missing: %s' % sorted(missing))
                sys.stderr.write('Wrong size: %s' % sorted(wrongsize))
            else:
                missing_tmp = set(phenos.keys()).difference(bcps)
                wrongsize_tmp = [bcp for bcp in bcps if actual_sizes[bcp] != expected_sizes[bcp]]
                missing = set([phenos[k] for k in missing_tmp])
                wrongsize = set([phenos[bcp] for bcp in wrongsize_tmp])
                if check is False:
                    wrongsize_tmp = []
                    wrongsize = set()
                sys.stderr.write('Missing: %s' % sorted(missing_tmp))
                sys.stderr.write('Wrong size: %s' % sorted(wrongsize_tmp))

            missing = ' '.join(map(str, [x + 1 for x in sorted(map(int, missing.union(wrongsize)))]))

            print('Re-running necessary tasks')
            rerun_cmd = os.path.join(prog_path, 'fastlmm_pipeline.py -u -s %s -m %s %s %s --tasks %s' % (species, memory, folder_base, ['', '-i'][igv], missing))
            print(rerun_cmd)
            sys.stdout.flush()
            if not dryrun:
                subprocess.call(rerun_cmd, shell=True, executable='/bin/bash')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--numeric', dest='numeric', help='expect numeric phenotype names',
                        default=False, action='store_true')
    parser.add_argument('-m', '--memory', dest='memory', help='megabytes of memory requested per job',
                        default=1200, action='store', type=int)
    parser.add_argument('-c', '--check-file-sizes', dest='check', help='check file sizes of binary output files when deciding if job is finished running',
                        default=False, action='store_true')
    parser.add_argument('-i', '--igv_output', dest='igv', help='if IGV output was selected, check for presence of .gwas files instead of .bcp files',
                        default=False, action='store_true')
    parser.add_argument('-s', '--species', dest='species', help='mouse or human',
                        default='mouse', action='store', choices=['human', 'mouse', 'dog', 'horse', 'cow', 'sheep'])
    parser.add_argument('-d', '--dryrun', dest='dryrun', help='dry run',
                        default=False, action='store_true')
    parser.add_argument('folders', nargs='+', help='results folder(s) to merge')
    args = parser.parse_args()
    folders = args.folders
    numeric = args.numeric
    memory = args.memory
    check = args.check
    igv = args.igv
    species = args.species
    dryrun = args.dryrun

    for folder in folders:
        if os.path.exists(folder):
            process(folder, numeric, memory, check, igv, species, dryrun)
        else:
            sys.stderr.write('Path not found: %s' % folder)
