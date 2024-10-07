# longRead Assembly pipeline 

**README for Nephele2 HPC USERs** \
Angelina Angelova \
Mar 2024 \



## List of abbreviations in text

* **SampName**: This stands for _Sample Name_ (not to be confused with _File Name_). Sample Names should be provided in the _$SampList_ and should _not_ contain file name extensions or tags (e.g. L001, \_R1, \_R2, .fastq.gz)
* **FileName**: Names of the files for the dataset. These contain tags and name extensions such as  L001, \_R1, \_R2, .fastq.gz, etc. 
* **ScriptName**: The name of the script wished to be executed (and they are numbered for convenience)
* **OutDir**: the output dir that will be created by script, containing the output files
* **suppPATH**: this is a parameter within the scripts that gives the path to supplementary scripts and data files used within the script. This path needs to be updated within each script prior to its **first use** (only) for each the machine/user that will be executing it
* **AA**: Amino acid
* **NT**: Nucleotide
* **bp**: base pairs
* **PWY**: pathway / pathway inference
* **TAX**: taxonomy / taxonomic




## Initial comments: 

This is a snakefile that is set to run & tested on Biowulf in a cluster setup => there are some rules, options and settings that are set up for biowulf HPC, that would need to be adjusted/amended for other HPCs or EC2.

1) Every tool that is used by the snakefile on Biowulf runs, is loaded as a module through the `cluster_config-biowulf.yaml`. For different HPC is to be used, create a different cluster_config file specific for that HPC.

2) The supplemental & R scripts used by the snakemake are provided under the supplemental/ folders of the pipeline, but some files need to be un-archived. Paths used in snakefile for each of those scripts & files, are set up in `config.yaml`, and need to be modified for the HPC in use. Note: the minpath tool provided with LoRA, is slightly modified from original. It is better to use the provided modified version

3) For the biowulf version of the pipeline, there was no global options settings. CLI user just manually unhash the block of paths associated with each classification DB (and comments out the un-needed DB paths).

4) The input fileNames of the `.fastq.gz` files for this snakemake read in through the mapping file. => The input fileNames read in follow a structure: `readsDIR + SampName + readSFFX + readsEXT`.


## Snakefile & supporting scripts

Scripts can be found on web at: [https://github.com/angelovaag/LoRA_public](https://github.com/niaid/microbiome-script-share/tree/main/lRA_pipeline) 

File structure is:

##### Main pipeline scripts: `snakefile-clust` , `config.yaml` && `utils.smk` 
##### Cluster setup scripts: `utils-clust.smk` & `submit_snake_biowulf.sh` 
#### Supporting, customized & R scripts: `supplementary`, 
#### Databases:  
* taxonomic & CARD databases: `TAXdb` files & `CARDdb` files can be downloaded using script \
`LoRA_public/supplementary/conda_installs/download_DB.sh`
* custom decontamination DBs: USER can either prepare their own custom `HOSTdb` files [using Kraken](https://telatin.github.io/microbiome-bioinformatics/Build-a-kraken2-database/) (add the —no-masking flag), or contact LoRA team for copy of the already built `HOSTdb` files linked below.

## Software tools, scripts, databases (from web)

#### assembly & assembly QC tools
* [Flye](https://github.com/fenderglass/Flye) (v3.15.3) ## an assembly tool
* [minimap2](https://github.com/lh3/minimap2)(v2.26) ## long-read alignment 
* [bbtools](https://jgi.doe.gov/data-and-tools/bbtools/) (v39.06) ## read alignment stats
* [samtools](http://www.htslib.org/) (v1.13) ## read alignment processing

#### gene & pathway analysis tools
* [Prodigal](https://github.com/hyattpd/Prodigal) (v2.6.3)
* [RGI](https://github.com/arpcard/rgi) (v6.0.2) ## Resistance Gene Finder
* [seqtk](https://github.com/lh3/seqtk) (v1.4)  ## seq extraction tool
* [eggNOG-mapper2](https://github.com/eggnogdb/eggnog-mapper) (v2.1.2)
* [subread](https://subread.sourceforge.net/) (v2.0.3)
* [MinPath](https://github.com/mgtools/MinPath) (v1.5) ## install, but use custom. Ask AA for details!! 

#### taxonomic classification tools
* [kraken2 && backen](https://github.com/DerrickWood/kraken2) (v2.1.2)   ## tax classification & decontamination tool
* [kronaTools](https://github.com/marbl/Krona/tree/master/KronaTools) (v2.8.1)
* [GTdb-tk](https://ecogenomics.github.io/GTDBTk/index.html) (v2.3.2) ## for binTAX classifications (might require singularity)

#### binning & binQC tools
* [CheckM](https://github.com/Ecogenomics/CheckM) (v1.1.3) ## bin QC
* [metabat](https://bitbucket.org/berkeleylab/metabat/src/master/) (v2.15)
* [blobTools](https://github.com/DRL/blobtools) (v1.1.1) ## for binTAX visualizations 

#### statistical and other tools 
* [pigz](https://zlib.net/pigz/) (v2.4)
* [R language](https://www.r-project.org/) (v4.3)
* [biom-format](http://biom-format.org/) (v2.1.6) ## making biom files






## Input requirements & files, at start

1) **mapping file**, set through `map_file` variable in `config.yaml` file. 

* The mapping file's fist column is used by the snakefile to read-in the Sample Names for initiation of pipeline & processing.
* The mapping file's second column is read in, but I couldnt make snakemake use those file names because it couldnt match them to sample names => it is not used
* The mapping file's third column is used  in the `DivPlots` steps to obtain groupings for the plotting steps (end steps)
* Mapping file should be of structure:

|#SampleID | FastqFileName               | TreatmentGroup |
|----------|-----------------------------| ---------------|
|Sample1   | Sample1_someSuffix.fastq.gz | Group1         |
|Sample2   | Sample2_someSuffix.fastq.gz | Group2         |
|Sample3   | Sample3_someSuffix.fastq.gz | Group2         |
|Sample4   | Sample4_someSuffix.fastq.gz | Group1         |

2) **reads DIR **: set through `readsDIR` in `config.yaml` to describe file names. Snakemake reconstructs file name structure with:  `readsDIR` + `SampName` + `readSFFX` + `readsEXT`. Reads could be:

3) **reads**: can be `fastq.gz` or `fastq` (theoretically, not tested). Reads can also be:

* _RAW reads_:  non-QCd reads, raw reads (still its better to have primers and adaptors removed from the sequence but sequences do not need to be HQ); or 
* _HQ20reads_: reads that have passed through minionQC pipeline (and are thus with quality >20). Those may have a `readSFFX` like "\_trimmed", "\_filtered", "\_pc_trimmed" or anything). `readSFFX` is a variable of a string that is not wished propagated as part of the SampName in output files; or
* _CORRreads_: reads that have passed not only QC but also a error-correction tool such as Canu or something (some users might have)


## Snakemake files and needed adjustments for initial run

At first use of the snakemake, adjust the following files (for the HPC or project settings needed)

1) **snakefile**: `snakefile-clust`. Adjust line 8 to point to `utils-clust.smk` (it can be generalized utils-clust.smk location)
	
2) **utils file**: `utils-clust.smk`. Already set for biowulf. Only need to adjust the SLURM variables, if switching HPCs

3) **cluster config file**: `cluster_config-biowulf.yaml`. Set up for Biowulf, need to change settings for each rule if switching HPCs. Also, if you dont like the log or job _name_ of a rule, you can change it here.

4) **project config file**: `config.yaml`. Adjust the `scritsDIR` & `clusterFILE` variables to correct paths. There is also `permanent PATHS` section that need to be set for each user. For project-related settings see **Pipeline USER options** section

5) **submit file**: `submit_snake_biowulf.sh` This is your main snakemake submission script. Adjust variables `sciptsDIR`, `logDIR` (logDIR is set to name slogs/ but you can change it). Also feel free to change `--jobs 4` flag in the snakemake line, to as many jobs at a time you want submitted. This script would also need further adjustments if switching HPCs. Also,  I have only tested it running in interactive session, but theoretically it should also do sbatch session. Keep in mind that if the session meets time wall, the submit script cannot submit cluster jobs anymore.




# Pipeline options (config.yaml)

Paths to **scriptDIR** (path to folder where all the scripts are) && **clusterFILE** (full path and name of the cluster_config file) need to be present and set. Also, set the user-dependent **permanent PATHS** at end of file (or set to use global ones from the shared dir)

### Input SETTINGS

I have this section prepped in 2 sets (one fo ONT data, one for PacBio data), but each project will only have one section unhashed. The other is just a shortcut for testing & example. 

**MAIN_ASMB_DIR_NAME**: set as whatever you want the assembly folder name to be. Technically, it should be just "assemblies/", but during testing is convenient to have variable name.

**map_file**: mapping file path & name. SampleID & TreatmentGroup are used by snakemake.

|#SampleID | FastqFileName               | TreatmentGroup |
|----------|-----------------------------| ---------------|
|Sample1   | Sample1_someSuffix.fastq.gz | Group1         |
|Sample2   | Sample2_someSuffix.fastq.gz | Group2         |
|Sample3   | Sample3_someSuffix.fastq.gz | Group2         |
|Sample4   | Sample4_someSuffix.fastq.gz | Group1         |

**readsDIR**: path and DIR where reads live

**readsSFFX**: any extra string in the file names, beyond the SampName, that you do not wish propagated through the pipeline names (e.g. \_001 or \_trimmed )

**readsExt**: if the file name extension is `.fastq.gz` or just `fastq` or `fq` or `fq.gz`. Technically I have only tested on .fastq.gz but I theoretically fastq files should work too

**SEQ TYPE**: whether you have ONT or PC data. String has to match `ONT` (Oxford NanoPore Tech) or `PB` (PacBio)

**DATA_QUAL**: Assembler & read mapper can take in RAW (no QC, but preferably at least adaptors removed), HQ (qc'd) or corrected (CORR) data (e.g. error corrected with something like Canu)
 	
 	
### Project SETTINGS:

**hostCONF**: the confidence level for kraken during HOST decontamination

**taxCOMF**: the confidence level of kraken during TAX classifications

<!--**readTAX**: boolean [def: F]: read-based taxonomic classification and read scoring by Kraken. It is rather _abundance inaccurate_ for long reads, but it is good for quick taxonomic screening of a dataset.-->

**asmb_polish_steps**: integer [def: 2]. The iteration of polishing steps for Flye to perform, during assembly. Technically, the tool's default is 1. But I like polished assemblies, so lRA's default is 2. You can set less or more depending on time and preference (they dont actually take too long)

**PREDICT_FUNCTION**: boolean [def: T]. Sets the pipeline to continue to gene finding, annotation and path inference. In some cases (if you ran pipeline before), this is not needed => I have made this convenience OFF function

**annotation_type**: selects which annotation type you want used for pathway inference: _KO against KEGG_ or _EC against MetaCyc_

**runRGI**: boolean [def: F ]. Choose if you want resistance gene finder run



### MAGs SETTINGs
These settings depend on bin_wCheckM being set to TRUE.


**bin_w_checkM**: boolean [def: F], if you wish MAGs created, TAX assessed & QCd (with metabat2, checkm2 & checkm1) (checM1 based taxonomic classification is provided here)

**GTDBtk**: boolean [def: F]. If you wish taxonomy of the MAGs established with GTdb-tk (better than checkm1 tax classification)

**checkM_plots**: quality plots assessed with checkM1 (since checkm2 is implemented, only some plots are produced)

**blobPlots**: boolean [def: T]: only the best quality plots for MAGs, ever!

### Databases SETTINGS 
 
**HOSTdb (decontamination databases)**: unhash **one line** of following options:
	
	 - defaultDB  (human+mouse\_db) [default]
	 - mosquitoDB (mosq\_db)
	 - nematodeDB (nematode\_db)
	 - marineDB   (marine\_db)
	 - plantsDB   (plants\_db)
	 - primatesDB (primates\_db)
	 - rodentsDB  (rodents\_db)
Note: those are provided upon request or custom build by USER

**TAXclassificaton DBs**: unhash elect **one block** of the following options

	- referenceDB (RefSeqDB) [default]
	- mouseGutDB  (MGBCdb)   
	- genomeTaxDB (GTdb)
	


# Submit a job

I usually start a small interactive session (e.g. 4CPUs/16G MEM), but theoretically submitting as sbatch is also possible (not tested). I have set up the `submit_snake_biowulf.sh` script to intake \${1} as a flag for dry run. Thus, 

- if you want to do a dry run of the snake, do `./submit_snake_biowulf.sh -np` inside an interactive session. (Technically any nonsense in \${1}  position will be interpreted as 'dry run')
- if you want to do an actual run of the snake: `./submit_snake_biowulf.sh` is sufficient

The submit script will created a `logsDIR` in the project DIR, to output each cluster's job log. Each log will carry the name of the sample that it works on. Jobs on global rules will carry the name of the function of the rule. 




## Overview of pipeline workflow

### Step 1) Decontamination

* Tools: **kraken2**, with `--gzip-compressed` flag (but I think kraken will ignore it if data is not gz)
* Rules: `decontam`
* Resources: Low demand on memory, depending on DB chosen by user (`$HOST DB`): ~2x size of DB. 
* Notes: 
  * To be performed against the `$HOST DB` chosen by user. 
  * for Biowulf, kraken is hardwired to expect zipped files. Flag `--gzip-compressed` may have to be removed if input is `.fastq`
  * for Biowulf, the input files are read in as: `readsDIR + sampName + readSFFX + readEXT`
  * Functionality: The script reads in `fastq.gz` files from whatever `readsDIR` is chosen and decontaminates against host DB chosen.
* **Expected Input**: 

  1) `$SampName` - the wildcard for each sample name 
  2) `$readsDIR` - input DIR chosen
  3) `$readSFFX` - suffix variable for the string used beyond the \$SampName & before the file extension ".fastq.gz" 
  4) `$readsEXT` - file name extension like `fastq.gz` 
* **Output**:
  * `decontam/SampName_{decontam.fastq.gz,_contamLOG.txt,_contamREPORT.txt}`: DIR with decontaminated seq ,log & report files (3x sample)

### Step 2) Assembly

* Tool: **Flye** with its `--nano-raw` or `--pacbio-raw` & `--meta` flags
* Rules: `assembly`
* Resources: High demand for time, memory & tmpDIR space. Usage builds up over time
* Functionality: assembles long reads
* Notes:
  * Creates a lot of large temp files & folders for each `$SampName`, so it is set to run in a `tmpDIR`
  * Only the final files (prefixed with `assembly`) are copied to the `$outDIR`
  * The `$outDIR` should be `assemblies/SampName_asmb/` by default, but ...
  * for Biowulf (SMEs & testing), the `$outDIR` main folder name is set with the `MAIN_ASMB_DIR_NAME` variable in `confif.yaml` file.
  * for Nephele, variable `MAIN_ASMB_DIR_NAME` should be set to **"assemblies/"** and _not be user option_
  * `$outDIR` will contain individual assembly subfolders for each sample, named `SampName_asmb/`
  * within each `$SampName_asmb/` folder, output file prefixes start with `assembly_`
* **Output**:
  * `assemblies/$SampName_asmb/assembly_graph.{gfa,gv}` (2 graph files)
  * `assemblies/$SampName_asmb/assembly{_info.txt,.fasta}` (2 main output files) 
  * `assemblies/$SampName_asmb/assembly.done,flye.log`  (2 log files)

### Step 3) Mapping decontaminated reads to assembly

**Info & Notes**:

* Tool: `minimap2`, `samtools`
* Resources: average
* Rules: `minimap_index`, `minimap2bam`
* Functionality: maps decontaminated reads to the assembled scaffolds. Then uses the mapped reads file (sam/bam) to enumerate reads mapped for each scaffold
* **Output**:
  * output within the `assemblies/$SampName_asmb/` folder
  * `assemblies/$SampName_asmb/assembly.{bam,bam.bai,mapLOG}`
 

### Step 4) Getting other assembly stats

**Info & Notes**:

* Tools: **samtools**, **metabat**, **bbtools**
* Rules: `getDEPTHs`, `asmbSUMMARY`
* Resources: average
* Functionality: Uses the assembly FASTA, BAM & INFO files, to calculate and add more stat columns to the assembly SUMMARY file from each sample. These stats are used later in the workflow
* **Output**:
  * `assemblies/$SampName_asmb/assemblySUMMARY.txt`



### Step 5) Taxonomic classifications

* Tools: **kraken2**, **kronatools**, **taxonkit**
* Rules: `asmbTAX_Kracken2`, `asmbTAX_addTAX_to_sumFILE`, `asmbTAX_kronaText`, `asmbTAX_kronaPlot`
* Resources: High memory depending on `$classDB` chosen 
* Functionality: assigns taxonomy to assembled scaffolds & uses assembly stats file (from step4) to estimate abundance. Also creates .html plots 
* **Output**:
  * `assemblies/$SampName_asmb/asmbTAX/$TAXdbNAME/{classLOG,$SampName{_4krona,_taxREPORT}}.txt`
  * `assemblies/$SampName_asmb/assemblySUMMARY_$dbNAME.txt`
  * `profilesTAX/$TAXdbNAME/kronabin/$SampName_4krona.txt` (Copies of _4krona files from each sample)
  * `profilesTAX/$TAXdbNAME/asmbTAX_plots.html` **---> for Nephele RESULTS page** 



### Step 6) Taxonomic profiling & visualizations

* Tools: **custom R and python3 TAX merging & plotting scripts**, **biom-format**,  
* Rules: `asmbTAX_mergedTAB`, `asmbTAX_biomTABs` & `asmbTAX_divPLOTs` (if len(samples)>1)
* Resources: low
* Functionality: takes in files from `profilesTAX/$TAXdbNAME/kronabin` folder, merges them into a community matrix & creates diversity plots. The merging script is set to if the R-version errors out, the python version is automatically attempted to produce  the same resulting table. 
* **Additional input**: 
  * the mapping file set with `$map_file` variable in `config.yaml` file
* **Output**:
  * `profilesTAX/$TAXdbNAME/merged_asmbTAX_{table_wLIN.txt,_json.biom}`
  * `profilesTAX/$TAXdbNAME/DivPlots/TAX_*.{pdf,png,txt,done}` (10 files)
  * `profilesTAX/$TAXdbNAME/DivPlots/TAX_BetaDiv_PCoA.png` **---> for Nephele RESULTS page** 



#### ----------- CONDITIONAL RULES ------------------------------

#### ----------- If `PREDICT_FUNCTION` == TRUE ------------------ [default = TRUE]
### Step 7) Gene finding & abundance scoring

* Tools: **prodigal**, **subread::featureCounts**, **eggnog-mapper**, **custom awk scripts**
* Rules: `asmbPWY_featurePrediction`, `asmbPWY_featureCounts`, `asmbPWY_featureTPM`, `asmbPWY_featureAnnotation`, `asmbPWY_geneSUMMARY`
* Resources: High demand on time & memory. 
* Functionality: predicts gene features, annotates genes and enumerates reads aligning to each gene/feature (basic gene stats)
* **Output**
  * `assemblies/$SampName_asmb/asmbPWY/genesSUMMARY.txt`
  * `assemblies/$SampName_asmb/asmbPWY/features/{feature.{faa,fna,gff,gtf},feature_{stats,counts_wTPM}.txt}`
  * `assemblies/$SampName_asmb/asmbPWY/annotations/annotations.emapper.{annotations,hits,seed_orthologs}`,
  * `assemblies/$SampName_asmb/asmbPWY/annotations.txt`

### Step 8) Gene extractions & gene stats gene table merging

* Tools: **custom AWK scripts**, **seqtk**
* Rules: `asmbPWY_geneSEQextractHDRS`, `asmbPWY_geneSEQextract`, `asmbPWY_aveGEN`, `asmbPWY_merge_aveGEN_files`
* Resources: low
* Functionality: uses information from gene finding step to extract gene sequences, basic gene stats and calculate TPM per gene and average per annotation.
* **Output**:
  * `assemblies/$SampName_asmb/asmbPWY/genes{EC,KO}.{faa,fna}`
  * `assemblies/$SampName_asmb/asmbPWY/$SampName_aveGEN.{ec,ko}.txt`
  * `profilesPWY/geneBIN{ec,ko}/$SampName_aveGEN.{ec,ko}.txt`
	
### Step 9) Inferring PWYs from gene information

* Tools: **minPath-PRNi**, **custom python3 scripts**
* Rules: `asmbPWY_prep_minPATH_DIR`, `asmbPWY_prep_minPATH_files`, `asmbPWY_run_minPath`,`asmbPWY_extract_minPATH_pathways`
* Resources: low
* Notes: minPath tool needs to be installed, but essentially our custom `minPATH-PRNi` script is used. It enables a manual set for MinPath's temporary files (with `--mps` flag), so these dont get overprinted from sample to sample. _thank you Poorani!_
* Functionality: use `annots_type` gene information (user election), to infer functional information from the genetic information. The steps use a few custom database files that map genes to pathways (PWYs) as well as the sample-specific gene and gene abundance information files. These files get temporarily copied into each sample's folder and be used by minPath locally, in order to avoid interference during processing. This means a few temp files are produced during these steps, but they are _not_ to be processed in `tmpDIR`
* **Output**: 
  * `assemblies/$SampName_amb/asmbPWY/pathways_$anno2pwy/minPATH.{log,report,details}.txt`  
  * `assemblies/$SampName_amb/asmbPWY/pathways_$anno2pwy/complete_pathways.txt`
	
### Step 10) Visualizing PWY profiles & merging PWY tables

* Tools: **kronatools**, **biom-format**, **custom R scripts**
* Rules:  `asmbPWY_run_genes2krona`, `asmbPWY_copy_4krona_text`, `asmbPWY_make_krona_pwyPLOT`, `asmbPWY_merge_PWYtables`, `asmbPWY_make_biom_PWYtables` ,`asmbPWY_divPLOTs`
* Resources: low
* Functionality: inferred PWY files are reformatted into visualization files & plotted. Tables are merged & diversity plots made for functional inference
* **Additional input**: 
  * the mapping file set with `$map_file` variable in `config.yaml` file
* **Outputs**
  * `assemblies/$SampName_amb/asmbPWY/pathways_$anno2pwy/$SampName_4krona.txt`
  * `profilesPWY/pwyBIN_$pwyDB/$SampName.txt` (copy of the `_4krona` file from each sample, no "\_4krona" string)
  * `profilesPWY/asmbPWYs_$pwyDB_plot.html` **---> for Nephele RESULTS page** 
  * `profilesPWY/merged_avePWY.{KEGG,MetaCyc}_{table,json}.{txt,biom}`
  * `profilesPWY/DivPlots_$pwyDB/PWY_*.{pdf,done)` (4 files)
	
	
### --------------- if `runRGI` == TRUE ----------------- [default = FALSE ]
### Step 11) Resistance Gene Identification

* Tools: **RGI** & **CARD_DB**
* Rules: `run_RGI`
* Resources: average
* Functionality: Uses features predicted in Step 7 & the `$CARD_db` to identify resistance genes
* **Output**:
  * `assemblies/$SampName_asmb/asmbPWY/RGI/rgiFNA_{main,raw}.{txt,json}` (3 files)
	
	
### ------------- if `bin_w_checkM` == TRUE ------------- [default = FALSE]
### Step 5.1) Binning of assembled scaffolds

* Tools: **metabat**
* Rules: `binning`
* Resources: average
* Functionality: scans scaffolds for specific features & separates them into separate multi-fasta files (called 'bins') that are assumed to be different organisms.
* **Outputs**:
  * `assemblies/$SampName_asmb/MAGs/mag.{[0-9]+,tooShort,unbinned}.fa` (number of files depends on sample)
  * `assemblies/$SampName_asmb/MAGs/{binsList.txt,binning.done}`

### Step 5.2) binQC workflow (CheckM)

* Tools: **m-tools** / **checkm**
* Rules: `checkm_lineage_wf`, `checkm_qa`, `checkm_coverage`, `checkm_profile`, 
* Resources: High. Very time & memory consuming steps. We hate CheckM
* Functionality: using checkM internal database, it determines the most likely taxonomic affiliations of the scaffolds within each bin, the level of completion of the most likely organismal affiliation of the entire bin, contamination level of the bin and other stats. 
* **Outputs**: 
  * `assemblies/$SampName_asmb/MAGs/checkM/{lineage.ms,checkm.log.lineageWF.done}`
  * `assemblies/$SampName_asmb/MAGs/mag{QA,TAX,COVR,PROFILES}.txt`

### Step 5.3) binQC summary & visualizations

* Tools: **kronatools**, **checkm**
* Rules: `join_profiles_tax`, `checkm_4krona`, `checkm_kronaPlot`, ~~`CheckM_plots (SME option only)`~~
* Resources: low
* Functionality: some processing of the CheckM outputs takes place for summary and visualization files.
* **Outputs**:
  * `assemblies/$SampName_asmb/MAGs/magSUMMARY.txt` 
  * `assemblies/$SampName_asmb/assemblySUMMARY_$dbNAME_wBINqc.txt` 
  * `assemblies/$SampName_asmb/MAGs/magTAXplots.html`
	


### ------------- if `GTDBtk` == TRUE -------------- [default = FALSE]
### Step 5.4) run GTDBtk	& produce summaries & visualizations

* Tools: **gtdb-tk**, **kronatools**,
* Rules: `binTAX_GTdbtk`, `binTAX_GTdb_4krona`, `binTAX_GTdb_kronaPlot`, `binTAX_add_GTdbTAX_to_covFILE`
* Resources: High! Time & memory & space demanding
* Functionality: in addition to checkM (completed in steps 5.{2,3}), GTDB-tk is run, to improve taxonomic resolution into TAX assignments of individual bins (the tax predictions are literally better!). Tool still only valid for prokaryotic organisms. The GTDB predictions are visualized & summarized.
* **Outputs**:
  * `assemblies/$SampName_asmb/MAGs/checkM/GTDBtk/` main folder with 4 sub-folders (`{align/,classify/,identify/,logs/}`)
  * `assemblies/$SampName_asmb/MAGs/checkM/GTDBtk/GTdb{-TAX.summary.tsv,-tk.done,_magTAXplots.html}` (3 summary fiiles)
	
### ----------- if `blobPlots` == TRUE ------------ [default = TRUE]
### Step 5.5) creating blob plots

* Tools: **blobtools**, **myblobtools script**
* Rules: `binTAX_blobtools`
* Resources: Time consuming but average on the computational resource consumption
* Functionality: We like this tool. It takes in the TAX classifications produced for each scaffold in step 5.0, the bins from step 5.1 and some stats from step 5.2, and creates nifty plots visualizing the content and quality of each MAG/bin. Tool does not depend on GTDBtk step
* **Outputs**: 
  * `assemblies/$SampName_asmb/MAGs/checkM/blobPlots/` main folder with 6x files per bin (lots of files).
  * `assemblies/$SampName_asmb/MAGs/checkM/blobs.done` (finished successfully)  




# end of file
