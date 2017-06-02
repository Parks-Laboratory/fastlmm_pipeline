# FastLMM Pipeline
[Documentation](http://microsoftgenomics.github.io/FaST-LMM/)

## Running the FASTLMM Pipeline on UW-Madison Clusters
This pipeline runs with the same inputs as [the epistasis pipeline](https://github.com/Parks-Laboratory/epistasis_pipeline), so 

## Preparing files for transfer, submit jobs (fastlmm_submit.py)
1. scp **data/**, **scripts/**, **fastlmm_submit.py** to submit server
	* **data/** contains
		**_prefix_.FILTERED.bim**,
		**_prefix_.FILTERED.bed**,
		**_prefix_.FILTERED.fam**,
		**_prefix_.FULL.bim**,
		**_prefix_.FULL.bed**,
		**_prefix_.FULL.fam**, and
		**_prefix_.pheno.txt**
	* **scripts/** contains
		**python.tar.gz** (portable Python 2.7 installation),
		**atlas.tar.gz** (ATLAS library), and 	
		**fastlmm_node.py**
1. `python fastlmm_submit.py <prefix> [options]`
	Use `ls results/<prefix> | wc -l` to check the current returned file numbers
1. run the following command on local machine: `scp -r <CONDOR_ADDRESS>:results <destination_directory_at_Parks_Lab>`

## Re-running jobs
The process is exactly the same as when [rerunning epistasis jobs](https://github.com/Parks-Laboratory/epistasis_pipeline/blob/redesign/master/README.md#re-running-jobs). 

Usage: Call **fastlmm_submit.py** with **--rerun _file_with_jobs_to_rerun_**, and whatever other flags and parameters are necessary to ensure these jobs succeed. Be careful not to change inputs/parameters in such a way that the failed-job number(s) wind up mapping to different work.
