#!/bin/bash
### Biowulf params
#SBATCH --job-name=wgsa3
#SBATCH --output=%x-%j.out
#SBATCH --partition=norm
#SBATCH --mail-type=ALL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --cpus-per-task=1            # CPUs
#SBATCH --mem=8G                     # Job memory request
#SBATCH --time=2-00:00:00
#SBATCH --export=NONE
#SBATCH --gres=lscratch:1
#SBATCH --signal=TERM@120

### This file is only slightly butchered for privicy. the gist is there

dryrun=${1}
if [[ ${dryrun} == "" ]]; then dryrun="" ; else dryrun="-np"; fi
echo "----> ${dryrun} <---- flag used (none if empty)"

#### I have not setup the cluster file for WGSA2 yet (thank you Katie for making it!), so I do general setup
# sbatchCMD="sbatch -c {cluster.threads} --mem={cluster.mem} \
# 			--output=${logDIR}/{cluster.log}-{cluster.jobname}--%j.txt \
# 			--partition={cluster.partition} --time={cluster.time} {cluster.extra}"
sbatchCMD='sbatch --cores 32 --mem=220g --time=24:00:00 --gres=lscratch:600 \
			--job-name={rule}-{wildcards.sampName}-{jobid} \
			--error=slogs/sLog_{rule}-{wildcards.sampName}-{jobid}.txt \
			--output=slogs/sLog_{rule}-{wildcards.sampName}-{jobid}.txt' ## -c {threads} --mem=196g
##{name} wildcard would be rule-name. Dont know how to call by clusterName 


runClusters=False ## [def: False]
if [[ ${runClusters} == "True" ]]; then 
	njobsFlag="4" ## 4 separate clusters will run based ont he sbatchCMD settings
	clustFlag=(--cluster "${sbatchCMD}" )
	echo "----> Cluster setting on. Running jobs in seprate clusters"
else 
	njobsFlag=${SLURM_CPUS_PER_TASK} ## everything will run in local interactive session using all resources
	clustFlag=()
	echo "----> nonCluster setting on. Running jobs in interactive session."
fi

useConda=False ##False ## [def: False]
condaPREFIX=/dbs/software/miniconda/4.10.3/envs/
if [[ ${useConda} == "True" ]]; then
    envFlag="--use-conda  --conda-prefix ${condaPREFIX}" ##" --conda-frontend mamba"
    echo "----> using CONDA environments"
else
  envFlag="--use-envmodules"
  echo "----> using modules"
  # conda activate ${condaPREFIX}/LoRA
fi

scriptsDIR=/MyScripts/minION_pipeline/longRead_snk-clust

wDIR=$(pwd)
## output wDIRectory; subdirectory for each sample
## will be created underneath where the output files will go
logDIR=${wDIR}/slogs 
mkdir -p ${logDIR}/

## project config file
configFILE=config.yaml    ### ${scriptsDIR}/config.yaml

### COMMENT OUT below the section that you DO NOT need ############
##### commands for biowulf  ##############
# some module installed to shared directory
# export MODULEPATH=/dbs/modulefiles:$MODULEPATH
module load snakemake ## /7.30.1
clusterFILE=${scriptsDIR}/cluster_config-biowulf.yaml

snakemake -s ${scriptsDIR}/snakefile-clust  -pr --rerun-triggers mtime \
	--cores $SLURM_CPUS_PER_TASK --configfile ${configFILE}  \
	--jobs ${njobsFlag} "${clustFlag[@]}" \
	--latency-wait 45   --max-jobs-per-second 1   \
	--nolock --keep-going --keep-incomplete ${envFlag} ${dryrun} 
