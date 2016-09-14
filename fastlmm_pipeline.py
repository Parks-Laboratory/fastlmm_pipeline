#!/usr/bin/env python

"""
Runs Plink, Generates Submit Files and Processes Inputs/Outputs to Cluster
"""

import sys
import os
import stat
import re
import subprocess
from math import log10, ceil
import textwrap
import operator
from datetime import datetime
import time
import re



debug = False

# input/output folders
root = os.path.split(os.path.realpath(sys.argv[0]))[0]
dataLoc = os.path.join(root, 'data')
condor_output_root = os.path.join(root, 'condor_out')
job_output_root = os.path.join(root, 'results')

# script locations
prog_path = os.path.join(root, 'scripts')
fastlmm_script = os.path.join(prog_path, 'fastlmm_wrapper.py')
merge_script = os.path.join(prog_path, 'merge_fastlmm_output.py')
plink_script = os.path.join(prog_path, 'plink')


# filter out SNPs with MAF < 5%, missing genotype frequency > 10%; convert to binary format
make_bed_cmd = '%(plink_location)s --tfile %(tfile_prefix)s --allow-no-sex --maf 0.05 --geno 0.1 --make-bed --out %(tfile_prefix)s %(plink_species)s'

# calculate MAF and missing genotype frequency
make_maf_cmd = '%(plink_location)s --bfile %(tfile_prefix)s --allow-no-sex --out %(tfile_prefix)s %(plink_species)s --freq'
make_missing_cmd = '%(plink_location)s --bfile %(tfile_prefix)s --allow-no-sex --out %(tfile_prefix)s %(plink_species)s  --missing'
clean_cmd = 'for f in %(tfile_prefix)s.{nosex,log}; do if [ -e $f ]; then rm $f; fi; done'



# set appropriate number of chromosomes (include X and Y)
species_chroms = {'human':24, 'mouse':21}

available_datasets = {}



class Tee(object):
    def __init__(self, filename):
        self.logfile = open(filename, 'a')

    def send_output(self, s):
        sys.stderr.write('%s\n' %s)
        self.logfile.write('%s\n' %s)

    def close(self):
        self.logfile.close()

def timestamp():
    return datetime.strftime(datetime.now(), '%Y-%m-%d_%H-%I-%S')

def make_intervals(tasks):
    # find contiguous intervals
    tasks_int = sorted(map(int, tasks))
    diffs = map(operator.sub, tasks_int[1:], tasks_int[:-1])

    intervals = []
    begin = end = min(tasks_int)

    for i, task in enumerate(tasks_int[1:]):
        if diffs[i] == 1:
            end = task
        else:
            intervals.append((begin, end))
            begin = end = task

    intervals.append((begin, end))
    return intervals

# function to make sure that fids and iids match across files
def check_fids_iids(prefix):
    def get_fids_iids(fn, skip=0):
        f = open(fn)
        for i in range(skip):
            f.readline()
        #fid_iid = [tuple(x.strip().split('\t')[:2]) for x in f.readlines()]
        fid_iid = [tuple(x.strip().split()[:2]) for x in f.readlines()]
        f.close()
        return fid_iid

    fid_iid_tfam, fid_iid_pheno = map(get_fids_iids, ['%s.tfam' % prefix, '%s.pheno.txt' % prefix], [0, 1])
    if fid_iid_tfam == fid_iid_pheno:
        return True
    elif set(fid_iid_tfam) == set(fid_iid_pheno):
        # FID/IID pairs can be in different order, but warn in case this is not intended
        log.send_output('FID/IID pairs in %s.pheno.txt are not in the same order as in %s.tfam' % (prefix, prefix))
        return True
    else:
        missing = set(fid_iid_tfam).difference(fid_iid_pheno)
        extra = set(fid_iid_pheno).difference(fid_iid_tfam)

        if missing:
            # OK to have more individuals in .tfam file than .pheno.txt file
            log.send_output('FID/IID pairs in %s.tfam but not in %s.pheno.txt:' % (prefix, prefix))
            log.send_output('\n'.join(['%s\t%s' % x for x in sorted(missing)]))
            return True

        if extra:
            log.send_output('FID/IID pairs in %s.pheno.txt but not in %s.tfam:' % (prefix, prefix))
            log.send_output('\n'.join(['%s\t%s' % x for x in sorted(extra)]))

    return False



# function to run plink
def populate_available(dataLoc, numeric, species):
    os.chdir(dataLoc)
    candidates = [x for x in os.listdir(dataLoc) if x.endswith(('.tped', '.tfam', '.pheno.txt', '.covar.txt'))]
    tped, tfam, pheno, covar = map(lambda suffix:set([x.replace(suffix, '') for x in candidates if x.endswith(suffix)]), ['.tped', '.tfam', '.pheno.txt', '.covar.txt'])
    geno = set(tped).intersection(set(tfam))
    plink_species = ['', '--%s' % species][species != 'human']

    for tfile_prefix in geno:
        sentinel = os.path.join(dataLoc, '%s.filtered' % tfile_prefix)
        tped = os.path.join(dataLoc, '%s.tped' % tfile_prefix)
        tfam = os.path.join(dataLoc, '%s.tfam' % tfile_prefix)
        local_numeric = numeric

        # check to see if genotypes have already been filtered
        if os.path.exists(sentinel) and \
                os.stat(sentinel)[stat.ST_MTIME] >= os.stat(tped)[stat.ST_MTIME] and \
                os.stat(sentinel)[stat.ST_MTIME] >= os.stat(tfam)[stat.ST_MTIME]:

            f = open('%s.bim' % tfile_prefix)
            n_snps = len(f.readlines())
            f.close()
        else:
            plink_location = plink_script
            for cmd_template in (make_bed_cmd, make_maf_cmd, make_missing_cmd, clean_cmd):
                cmd = cmd_template % locals()
                log.send_output(cmd)
                subprocess.check_call(cmd, shell=True, stderr=subprocess.STDOUT)
        
            chroms = map(str, range(1, species_chroms[species] + 1))
            n_snps = int( subprocess.Popen(['wc', '-l', os.path.join(dataLoc, '%s.bim' % tfile_prefix)], stdout=subprocess.PIPE).communicate()[0].split()[0] )

        sys.stdout.flush()

        if tfile_prefix in pheno:
            f = open('%s.pheno.txt' % tfile_prefix)
            # skip FID and IID columns
            headers = f.readline().strip().split('\t')[2:]
            pheno_data = f.readlines()
            f.close()

            # if FID/IID pairs match between .tfam and .pheno.txt, proceed
            # if not, go to next file
            if check_fids_iids(tfile_prefix):
                pass
            else:
                continue

            # check that phenotype names don't have spaces
            problematic = [x for x in headers if ' ' in x]

            # check that phenotypes exist and have unique names
            if not problematic and headers and len(set(headers)) == len(headers):
                # create phenotype key file for later reference
                f = open(os.path.join(dataLoc, 'pheno.key.txt'), 'w')
                fmt = '%%0%sd' % int(ceil(log10(len(headers))))
                f.write( '\n'.join(['%s\t%s' % (fmt % i, phen) for i, phen in enumerate(headers)]) )
                f.close()

                available_datasets[tfile_prefix] = {'n_pheno': len(headers),
                                                    'n_indivs': len(pheno_data),
                                                    'n_snps': n_snps}
                if tfile_prefix in covar:
                    available_datasets[tfile_prefix]['covar'] = '%s.covar.txt' % tfile_prefix
                else:
                    available_datasets[tfile_prefix]['covar'] = None

                # check numeric option
                # if headers are supposed to already be numeric, verify
                # if verification fails, fall back to regular numeric behavior
                if local_numeric == 2:
                    try:
                        numeric_headers = map(int, headers)
                    except ValueError:
                        log.send_output('Non-numeric values in phenotype names for %s; falling back to numeric=1 behavior' % tfile_prefix)
                        local_numeric = 1

                available_datasets[tfile_prefix]['numeric'] = ['', '-n %s' % local_numeric][numeric > 0]

            elif len(set(headers)) < len(headers):
                duplicates = sorted([x for x in set(headers) if headers.count(x) > 1])
                log.send_output('Duplicated phenotype names found in %s.pheno.txt:' % tfile_prefix)
                log.send_output('\n'.join(map(str, duplicates)))
            elif problematic:
                log.send_output('Spaces exist in the following phenotype names:')
                log.send_output('\n'.join(sorted(problematic)))
                log.send_output("cat <(head -1 %(dataLoc)s/%(tfile_prefix)s.pheno.txt | sed 's/ //g') <(tail -n+2 %(dataLoc)s/%(tfile_prefix)s.pheno.txt) > tmp.txt && mv -f tmp.txt %(dataLoc)s/%(tfile_prefix)s.pheno.txt" % locals())
            else:
                log.send_output('No phenotypes found in %s.pheno.txt!' % tfile_prefix)


# function to submit to clusters
def process_all(covar=False, memory=1024, species='mouse', maxthreads=1, featsel=False, exclude=False, condition=None):
    '''Processes all datasets'''
    process(sorted(available_datasets.keys()), covar=covar, memory=memory, species=species, maxthreads=maxthreads, featsel=featsel, exclude=exclude, condition=condition)

def process(datasets, covar=False, memory=1024, tasks=None, species='mouse', maxthreads=1, featsel=False, exclude=False, condition=None):

    os.chdir(root)

    # submit and executable files
    submit_template = textwrap.dedent(
    '''# FastLMM Submit File

    universe = vanilla
    log = %(condor_output)s/fastlmm_$(Cluster).log
    error = %(condor_output)s/fastlmm_$(Cluster)_$(Process).err

    InitialDir = %(root)s/results/%(dataset)s
    executable = %(root)s/fastlmm_%(dataset)s.sh
    arguments = $(Process)
    output = %(condor_output)s/fastlmm_$(Cluster)_$(Process).out

    should_transfer_files = YES
    when_to_transfer_output = ON_EXIT
    transfer_input_files = %(root)s/fastlmm_%(dataset)s.sh,%(prog_path)s/,%(dataLoc)s/

    request_cpus = 1
    request_memory = %(use_memory)sMB
    request_disk = 1GB

    queue %(n_pheno)s
    ''').replace('\t*', '')

    exec_template = textwrap.dedent(
    '''#!/bin/bash

    # untar your Python installation
    tar -xzvf python.tar.gz

    # make sure the script will use your Python installation
    export PATH=$(pwd)/python/bin:$PATH

    # run your script
    python fastlmm_wrapper.py %(covFile)s %(numeric)s %(igv)s %(debug)s %(species)s %(maxthreads)s %(feature_selection)s %(exclude)s %(condition)s %(dataset)s $1 >& fastlmm_wrapper.py.output.$1
    ''').replace('\t*', '')

    submission_cmd = 'condor_submit fast_lmm.sub'

    # generate submit files and submit to cluster
    for dataset in datasets:
        if dataset in available_datasets:
            params = available_datasets[dataset]

            # set memory and max threads
            if memory is None:
                # try to guess optimal amount of memory to request
                # regression coefficients determined empirically; fudge factor added
                local_memory = min(4096*16, int(1.65 * params['n_indivs'] * params['n_snps'] / 200000) + 440 + 300)
            else:
                local_memory = min(4096*16, memory)

            total_memory = local_memory
            if local_memory > 4096:
                maxthreads = max(maxthreads, min(16, int(ceil(local_memory / 4096.)) + 0))
                local_memory = 4096


            # generate output files
            condor_output = os.path.join(condor_output_root, dataset)
            if not os.path.exists(condor_output):
                os.makedirs(condor_output)

            job_output = os.path.join(job_output_root, dataset)
            if not os.path.exists(job_output):
                os.makedirs(job_output)

            if covar and params['covar'] is None:
                log.send_output('Specified --covar but no covariate file exists; ignored')

            if condition:
                condition = condition[0]

            params.update({'root': root,
                           'dataLoc': dataLoc,
                           'dataset': dataset,
                           'job_output': job_output,
                           'condor_output': condor_output,
                           'covFile': ['', '-c %s' % params['covar']][covar and params['covar'] is not None],
                           'fastlmm_script': fastlmm_script,
                           'debug': ['', '--debug'][debug],
                           'prog_path':prog_path,
                           'timestamp':datetime.ctime(datetime.now()),
                           'species': '-s %s' % species,
                           'maxthreads':'--maxthreads %s' % maxthreads,
                           'feature_selection':['', '--feature-selection'][featsel],
                           'exclude':['', '--exclude'][exclude],
                           'condition': ['', '--condition %s' % condition][condition is not None],
                           'use_memory': local_memory})

            maxthreads_option = ['', '-pe shared %s' % maxthreads][maxthreads > 1]

            if tasks is None:
                intervals = [(1, params['n_pheno'])]
            else:
                intervals = make_intervals(tasks)


            for low, high in intervals:

                submit_file = open( 'fastlmm_%(dataset)s.sub' % params, 'w')
                submit_file.write( (submit_template % params).replace(',,', ',') )
                submit_file.close()

                exec_file = open( 'fastlmm_%(dataset)s.sh' % params, 'w')
                exec_file.write( exec_template % params )
                exec_file.close()

                subprocess.call('chmod +x fastlmm_%(dataset)s.sh' % params, shell = True)

                # submit jobs to condor
                condor_cluster = subprocess.Popen(['condor_submit', 'fastlmm_%(dataset)s.sub' % params], stdout=subprocess.PIPE).communicate()[0]
                condor_cluster = re.search('\d{4,}', condor_cluster).group()
                print("Submitting Jobs to Cluster %s" % condor_cluster)
                log.send_output("%s was sent to cluster %s at %s" % (params['dataset'], condor_cluster, timestamp()))

                # check when jobs are done
                #n_running = 1
                #while n_running > 0:
                    #check_status = subprocess.Popen(['condor_q'], stdout=subprocess.PIPE).communicate()[0]
                    #subprocess.call('condor_q', shell = True)
                    #n_running = check_status.count(condor_cluster)
                    #time.sleep(120)

                #subprocess.call('mv *gwas %(job_output)s/' % params, shell = True)
                #subprocess.call('mv *output* %(condor_output)s' % params, shell = True)
                
                #subprocess.call('mv fastlmm*sub %(condor_output)s' % params, shell = True)
                #subprocess.call('mv fastlmm*sh %(condor_output)s' % params, shell = True)

                #log.send_output("%s finished at %s" % (params['dataset'], timestamp()))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='''Runs FaST-LMM on datasets found in specified location (looks in %s by default). Each dataset should have PLINK-formatted genotype (.tped, .tfam) and alternate phenotype (*.pheno.txt) files.  Optional covariate files should be named with the same prefix as the other files and end in .covar.txt .
''' % dataLoc)
    #parser.set_usage('''%(prog)s [options] [dataset1] [dataset2] ... (runs all datasets if unspecified)
    #PLINK-formatted genotype (*.tped, *.tfam) and alternate phenotype (*.pheno.txt) files should be placed in ''' + dataLoc)
    parser.add_argument('-l', '--list', dest='list_dataset', help='lists datasets to process, does not do processing',
                        default=False, action='store_true')
    parser.add_argument('-d', '--datadir', dest='datadir', help='specifies folder to search for raw data',
                        default=dataLoc, action='store')
    parser.add_argument('-o', '--outputdir', dest='outputdir', help='specifies output folder',
                        default=job_output_root, action='store')
    parser.add_argument('-c', '--covar', dest='covFile', help='use covariate file',
                        default=False, action='store_true')
    parser.add_argument('-s', '--species', dest='species', help='mouse or human',
                        default='mouse', action='store', choices=['human', 'mouse', 'dog', 'horse', 'cow', 'sheep'])
    parser.add_argument('-m', '--memory', dest='memory', help='amount of RAM (in megabytes) requested per job',
                        default=None, action='store', type=int)
    parser.add_argument('--maxthreads', dest='maxthreads', help='maximum # of threads to use',
                        default=1, action='store', choices=range(1, 17), type=int)
    parser.add_argument('-f', '--feature-selection', dest='featsel', help='perform feature selection',
                        default=False, action='store_true')
    parser.add_argument('-e', '--excludeByPosition', dest='exclude', help='exclude SNPs within 2Mb of tested SNP from kinship matrix construction',
                        default=False, action='store_true')
    parser.add_argument('-n', '--numeric_phenotype_id', dest='numeric', help='convert phenotype names to numbers (for safety)',
                        nargs='?', default=0, const=1, type=int, action='store', choices=[0, 1, 2])

    parser.add_argument('-q', '--quiet', dest='debug', help="suppress debugging output",
                        default=True, action='store_false')
    parser.add_argument('dataset', metavar='dataset', nargs='*', type=str, help='dataset(s) to process')
    parser.add_argument('--tasks', dest='tasks', metavar='TASK', nargs='+', help='run only specified sub-tasks (specify only one dataset when using this option)', type=int)
    parser.add_argument('--condition', dest='condition', help='condition on SNP {snp_id}',
                        action='store', nargs=1)

    args = parser.parse_args()
    dataLoc = args.datadir
    job_output_root = args.outputdir
    list_data = args.list_dataset
    covFile = args.covFile
    memory = args.memory
    datasets = args.dataset
    numeric = args.numeric 
    species = args.species.lower()
    maxthreads = args.maxthreads
    featsel = args.featsel
    exclude = args.exclude
    debug = args.debug
    tasks = args.tasks 
    condition = args.condition

    if debug:
        log = Tee('fastlmm_pipeline-%s.log' % timestamp())
    else:
        log = Tee('/dev/null')

    if tasks and len(datasets) > 1:
        log.send_output('More than one dataset specified along with --tasks option; quitting')
        log.close()
        sys.exit(0)

    log.send_output('Searching for raw data in %s' % dataLoc)
    populate_available(dataLoc, numeric, species)

    if list_data:
        log.send_output('\nAvailable datasets:')
        if available_datasets:
            log.send_output('\n'.join(sorted(available_datasets.keys())))
        else:
            log.send_output('\tNone')
    else:
        if datasets:
            log.send_output('\n'.join(sorted(datasets)))
            process(sorted(datasets), covFile, memory, tasks, species=species, featsel=featsel, exclude=exclude, condition=condition)
        else:
            log.send_output('\n'.join(sorted(available_datasets.keys())))
            process_all(covFile, memory, species=species, featsel=featsel, exclude=exclude, condition=condition)

    log.close()


