#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=16g
#SBATCH --job-name taxkit
#SBATCH --mail-type END
#SBATCH -o slogs/slog_%x_%j
#SBATCH --time=12:00:00
#SBATCH --partition=norm  
# if [ ! -d "slogs" ]; then mkdir slogs; fi

# module load taxonkit
## converts TAXid to  Lineage & TAXranks
## using kronaTAX files (1 or many samples [taxID SampName... ])

inFILE=${1%%.*}
# TAXDB=${2} ;  # manual use
# 			 if [[ ${TAXDB} == "" ]] ;            then TAXONKIT_DB=/usr/local/apps/taxonkit/TAXONKIT_DB ; fi # default from AUG 2022
# 			 if [[ ${TAXDB} == "plusPFPV" ]] ;    then TAXONKIT_DB=/data/BCBB_microbiome_db/blobtools-data/; fi ## JUL 2023
# 			 if [[ ${TAXDB} == "GTdb" ]] ;        then TAXONKIT_DB=/data/BCBB_microbiome_db/Kraken2_db/taxonomic_db/GTdb_r207/GTdb-taxdump_JUL2023/; fi
# 			 if [[ ${TAXDB} == "MBGCdb" ]];       then TAXONKIT_DB=/data/BCBB_microbiome_db/Kraken2_db/taxonomic_db/MGBCdb/taxonomy ; fi 
												   #TAXONKIT_DB=/data/BCBB_microbiome_db/Kraken2_db/db_prep_files/taxonomy/ ## from 2021; 
TAXONKIT_DB=${2} ## set up through snakemake
			
taxIDfield=${3} ; if [[ ${taxIDfield} == "" ]]; then taxIDfield=1 ; fi 
threads=${4} ; if [[ ${threads} == "" ]]; then threads=$SLURM_CPUS_PER_TASK ; fi
tmp_dir=${5} ; if [[ ${tmp_dir} == "" ]] ; then tmp_dir=tmp_dir/ ; fi

outFILE=${inFILE}_wLIN.txt


mode=kronaTAXfile
taxonkit lineage --taxid-field ${taxIDfield} ${inFILE}.txt -j ${threads} --delimiter ";" -o ${tmp_dir}/tmp_taxonkit1.txt --data-dir ${TAXONKIT_DB}
linField=$(awk -F'\t' '{print NF; exit}' ${tmp_dir}/tmp_taxonkit1.txt)
taxonkit reformat --lineage-field ${linField}   ${tmp_dir}/tmp_taxonkit1.txt --delimiter ";" -F \
	-f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" -R "Unclassified" -o ${tmp_dir}/tmp_taxonkit2.txt --data-dir ${TAXONKIT_DB}
cat ${tmp_dir}/tmp_taxonkit2.txt | sed 's/unclassified //g' | \
	sed 's/Candidatus //g' | sed 's/[][]//g' |  sed 's/cellular organisms.//g' | \
	sed '1 s|\t\t$|\tLineage\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies|'| \
	sed 's|\t\t|\tUnassigned\tUnassigned\tUnassigned\tUnassigned\tUnassigned\tUnassigned\tUnassigned\tUnassigned|g' | \
	sed 's|\tsuperkingdom\t.*$|\tUnassigned\tUnassigned\tUnassigned\tUnassigned\tUnassigned\tUnassigned\tUnassigned\tUnassigned|g' # >  ${outFILE} 


# rm -rf ${tmp_dir}

	## Note: the superkingdom sed: when taxonKit is missing a lineage, this line fills in some gaps
	# 	sed 's|superkingdom\t.*$|\tUnassigned\tUnassigned\tUnassigned\tUnassigned\tUnassigned\tUnassigned\tUnassigned\tUnassigned|g' 
	
#-r "Unassigned" will add "unnassigned" for each empty rank
#-R will include unclassified 

#great job this script does!!! but it works only based on taxID (NCBI tax)