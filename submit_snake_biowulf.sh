#!/bin/bash
### Biowulf params
#SBATCH --job-name=wgsa3
#SBATCH --output=%x-%j.out
#SBATCH --partition=norm
#SBATCH --mail-type=ALL         	 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --cpus-per-task=1            # CPUs
#SBATCH --mem=8G                     # Job memory request
#SBATCH --time=2-00:00:00
#SBATCH --export=NONE
#SBATCH --gres=lscratch:1
#SBATCH --signal=TERM@120

### Locus params
#$ -N wgsa3
#$ -cwd
#$ -j y
#$ -l h_vmem=8G
#$ -o .
#$ -m e
#$ -M your_address@niaid.nih.gov

dryrun=${1}
if [[ ${dryrun} == "" ]]; then dryrun="" ; else dryrun="-np"; fi
echo "----> ${dryrun} <---- flag used (none if empty)"


nonClust=True ## [def: False]
if [[ ${nonClust} == "True" ]]; then 
	clustFlag="--jobs ${SLURM_CPUS_PER_TASK}"
	echo "----> nonCluster setting on. Running jobs in interactive session."
else 
	clustFlag='--cluster "${sbatchCMD}" --jobs 4' 
	echo "----> Cluster setting on. Running jobs in seprate clusters"
fi

useConda=False ##False ## [def: False]
condaPREFIX=/data/BCBB_microbiome_db/software/miniconda/4.10.3/envs/
if [[ ${useConda} == "True" ]]; then
    envFlag="--use-conda  --conda-prefix ${condaPREFIX}" ##" --conda-frontend mamba"
    echo "----> using CONDA environments"
else
	envFlag="--use-envmodules"
	echo "----> using modules"
	# conda activate ${condaPREFIX}/LoRA
fi

scriptsDIR=/home/angelovaag/MyScripts/minION_pipeline/longRead_snk-clust

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
# export MODULEPATH=/data/BCBB_microbiome_db/modulefiles:$MODULEPATH
module load snakemake ## /7.30.1
clusterFILE=${scriptsDIR}/cluster_config-biowulf.yaml


sbatchCMD="sbatch -c {cluster.threads} --mem={cluster.mem} \
			--output=${logDIR}/{cluster.log}-{cluster.jobname}--%j.txt \
			--partition={cluster.partition} --time={cluster.time} {cluster.extra}"
##{name} wildcard would be rule-name. Dont know how to call by clusterName 


snakemake -s ${scriptsDIR}/snakefile-clust  \
	--jobname "{cluster.log}-{cluster.jobname}--{jobid}" \
	--configfile ${configFILE} ${clustFlag} --cluster-config ${clusterFILE} \
	--latency-wait 120 --max-jobs-per-second 1 --rerun-triggers mtime \
	--nolock --keep-going --keep-incomplete ${envFlag} ${dryrun}

## --jobname is the name on the bash dashboard
## --jobname "{cluster.jobname}-{jobid}" works, and is also %x
## --cluster "${sbatchCMD}" --cluster-config ${clusterFILE} --jobs 4 \

## --output is the name of the log file
## --output=${logDIR}/{cluster.log}__%x-%j.out.txt  would mean 
### ${logDIR}/{cluster.log}__{cluster.jobname}-{jobid}__{SLURM_JOB_ID} ==> too much


# ##### commands for locus  ##############
# export MODULEPATH=/hpcdata/bcbb/shared/microbiome_share/modules:$MODULEPATH
# module load snakemake/5.4.0-Python-3.6.7 || exit 1
# clusterFILE=${scriptsDIR}/fm.locus.yaml
# drmaacmd=" -l h_vmem={cluster.h_vmem} -j y -pe threaded {cluster.threads} {cluster.extra}"


# snakemake -s ${scriptsDIR}/functional_metagenomics.smk --jobs 32 --jobname "{name}.{cluster.jobname}.{jobid}" \
# 	  --configfile ${configFILE} \
# 	  --drmaa "${drmaacmd}" --cluster-config ${clusterFILE} --drmaa-log-dir ${outputdir}/{cluster.log} \
# 	  --nolock --keep-going --keep-incomplete --use-envmodules
