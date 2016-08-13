# FastLMM Pipeline

## Running the FASTLMM Pipeline on UW-Madison Clusters

1) add files (tfam, tped, pheno, covar) to the **fastlmm/data/** directory
2) copy **fixpheno.sh** file into the **fastlmm/data/** directory
3) run the command `sh fixpheno.sh PREFIX` where PREFIX is the actual prefix for the data set
4) move to the **fastlmm/ **directory
5) run the command `python fastlmm_pipeline.py --covar` (with additional arguments as necessary)
	Note: for a detailed list of additional arguments, run the command `python fastlmm_pipeline.py -h`
6) wait while pipeline runs
7) transfer the directory **results/** back to the server
8) check the **condor_out/** directory for any errors
9) clear out the **data/** directory for next time


NOTE:

* Sometimes the "condor_q" command fails and quits the program unexpectedly. In this case, rerun the pipeline
* Ensure that files scripts/fastlmmc, scripts/plink, fixpheno.sh are executable (green). If not run the command `chmod +x FILE_NAME`
