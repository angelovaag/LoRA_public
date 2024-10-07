###set conda main path
condaMAIN=${1}

### example: condaMAIN=data/BCBB_microbiome_db/software/miniconda/4.10.3/

### activate LoRA base
#source ${condaMAIN}/etc/profile.d/conda.sh && conda activate base
#LoRA_path=/data/BCBB_microbiome_db/software/miniconda/4.10.3/envs/LoRA/
# conda activate $LoRA_path

conda config --add channels defaults bioconda conda-forge biopython r ursky agbiome
conda install mamba

conda create -p ${condaMAIN}/envs/LoRA
conda create -p ${condaMAIN}/envs/RGI
conda create -p ${condaMAIN}/envs/checkm2
conda create -p ${condaMAIN}/envs/blobtools



mamba activate ${condaMAIN}/envs/LoRA
mamba install --prefix ${condaMAIN}/envs/LoRA/ pigz biopython kraken2 krakentools krona  \
	minimap2 metabat2 bbtools bbmap taxonkit subread prodigal eggnog-mapper seqtk biom-format \
	r-base r-essentials r-ggplot2 r-dplyr checkm-genome gtdbtk samtools flye boost=1.85.0
mamba update pigz biopython kraken2 krakentools krona  \
	minimap2 metabat2 bbtools bbmap taxonkit subread prodigal eggnog-mapper seqtk biom-format \
	r-base r-essentials r-ggplot2 r-dplyr checkm-genome gtdbtk samtools flye
# boost package is added for the proper libboost_program_options.so.1.85.0 install for metabat
mamba install --prefix ${condaMAIN}/envs/LoRA -y -c r \
r-fossil r-survival r-ape r-ampvis2  r-easypackages  r-docopt # r-limma  r-reshape 
## -y flag gives conda auto yes
mamba deactivate

mamba activate ${condaMAIN}/envs/rgi     && \
mamba install --prefix ${condaMAIN}/envs/rgi biopython rgi  && \
mamba update rgi && \
mamba deactivate


mamba activate ${condaMAIN}/envs/checkm2 && \
mamba install --prefix ${condaMAIN}/envs/checkm2 python biopython checkm2  && \
mamba update checkm2 && \
checkm2 database --download --path ${condaMAIN}/envs/checkm2
checkm2 database --setdblocation ${condaMAIN}/envs/checkm2/CheckM2_database/uniref100.KO.1.dmnd 
mamba deactivate


## add blobtools? with conda? 
mamba activate ${condaMAIN}/envs/blobtools
wget https://github.com/DRL/blobtools/archive/refs/tags/blobtools_v1.1.1.tar.gz -O ${condaMAIN}/envs/blobtools/blobtools_v1.1.1.tar.gz
tar -xzvf ${condaMAIN}/envs/blobtools/blobtools_v1.1.1.tar.gz -C ${condaMAIN}/envs/blobtools/
rm ${condaMAIN}/envs/blobtools/blobtools_v1.1.1.tar.gz
mamba install -y -c anaconda matplotlib docopt tqdm wget pyyaml git --prefix ${condaMAIN}/envs/blobtools
pip install pysam --user  ## mamba install --prefix ${condaMAIN}/envs/blobtools -c bioconda pysam  --update-deps

## installs in ~/.local/bin but needs env to run
conda deactivate

