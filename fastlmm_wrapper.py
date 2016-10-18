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
from fastlmm.util.runner import LocalInParts
import pysnptools.util
from pysnptools.util.pheno import loadOnePhen
from pysnptools.snpreader import Bed

root = os.path.split(os.path.realpath(sys.argv[0]))[0]
fastlmmc = './fastlmmc'

# ignore Y chromosome
species_chroms = {'human':23, 'mouse':20}


# run fastlmmc
def run_fastlmmc(dataset, output_dir, pheno_index, covFile=None, species='mouse', maxthreads=1, featsel=False, exclude=False, condition=None):
    
    # commands from fastlmmc:
    # maxthreads
    # condition
    # exclude by position
    
    # if condition:
    #     condition = '-SnpId1 %s' % condition[0]
    # else:
    #     condition = ''

    # temporary kludge because -excludeByPosition option is slow (at least for v2.05 and v2.06)

    bfile = dataset
    snp_reader = Bed(bfile)

    pheno_file = loadOnePhen('%s.pheno.txt' % dataset, i_pheno = pheno_index)
    phenotype_name = pheno_file['header'][0] 

    v = globals()
    chroms = map(str, range(1, species_chroms[species] + 1))
    v.update(locals())
    
    # loop through chromosomes and run
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
    v.update(locals())
    final.to_csv('%(output_dir)s/%(phenotype_name)s.gwas' % v, sep='\t', index=False)


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    #parser.set_usage('''%prog [options] dataset''')
    parser.add_argument('-s', '--species', dest='species', help='mouse or human',
                        default=None, action='store')
    parser.add_argument('-c', '--covariate_file', dest='covFile', help='use covariate file',
                        default=None, action='store')
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

    run_fastlmmc(dataset, output_dir, pheno_index, covFile, species, maxthreads, featsel=featsel, exclude=exclude, condition=condition)


