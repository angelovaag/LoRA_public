### set up a main path for all the DBs installations
DB_paths=${1}


## Decontam DBs download

### ------------

## RefSeq_plusPVFP db [~320Gb]:
RefSeq_DL_link=https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20240605.tar.gz
RefSeq_db=${DB_paths}/plusPVFP_Jun2024

mkdir -p ${RefSeq_db}
wget ${RefSeq_DL_link} -O ${RefSeq_db}/plusPFPV_kr2DB_Jun2024.tar.gz
tar -xzvf ${RefSeq_db}/plusPFPV_kr2DB_Jun2024.tar.gz.tar.gz -C ${RefSeq_db}/

## to get TAXDUMP, dont use taxonkit:  ## creates nonsense taxdumps
## taxonkit create-taxdump ${RefSeq_db}/ktaxonomy.tsv --out-dir ${RefSeq_db}/taxdump_PFPV_Jun2024 --threads 16
kraken2-build --download-taxonomy --db ${RefSeq_db}/taxdump_Jun2024 --threads 32

ln -s /usr/local/apps/kronatools/2.8.1/Krona-2.8.1/
./updateTaxonomy.sh ${RefSeq_db}/taxdump_PFPV_Jun2024/ --preserve --only-build
### ------------

## MGCB db [~40Gb]:
MGBC_db=${DB_paths}/MGBCdb_May2021/
MGBC_DL_link=https://zenodo.org/record/4836362/files/MGBC-26640_KrakenBracken.tar.gz?download=1

mkdir -p ${MGBC_db}
wget ${MGBC_DL_link} -O ${MGBC_db}/MGBCdb_May2021.tar.gz
tar -xzvf ${MGBC_db}/MGBCdb_May2021.tar.gz -C {MGBC_db}/

## to get TAXDUMP, dont use taxonkit:  ## creates nonsense taxdumps
## taxonkit create-taxdump ${MGBC_db}/ktaxonomy.tsv --out-dir ${MGBC_db}/taxdump_MGBCdb_May2021 --threads 16
## MGBCdb comes with its own taxdump, or you can download it from ??
ln -s ${MGBC_db}/taxonomy/ ${RefSeq_db}/taxdump_MGBCdb_May2021

## ln -s /usr/local/apps/kronatools/2.8.1/Krona-2.8.1/
./updateTaxonomy.sh ${RefSeq_db}/taxdump_MGBCdb_May2021 --preserve --only-build
### ------------

## CARD_db for RGI [~50Mb]
wget https://card.mcmaster.ca/latest/data -O ${DB_paths}/card-data.tar.bz2
tar -xzvf ${DB_paths}/card-data.tar.bz2 -C ${DB_paths}/CARD_db-v3.2.9/


## EGGNOG_DB for emapper: [~50G]
mkdir -p ${DB_paths}/eggnog_db_v5.0.2/

wget -p http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog.db.gz            -O ${DB_paths}/eggnog_db_v5.0.2/eggnog.db.gz
wget -p http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz      -O ${DB_paths}/eggnog_db_v5.0.2/eggnog.taxa.db.tar.gz
wget -p http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz -O ${DB_paths}/eggnog_db_v5.0.2/eggnog_proteins.dmnd.gz
tar -xzvf ${DB_paths}/eggnog_db_v5.0.2/eggnog.taxa.db.tar.gz   -C ${DB_paths}/eggnog_db_v5.0.2/eggnog.taxa.db
tar -xzvf ${DB_paths}/eggnog_db_v5.0.2/eggnog.db.gz            -C ${DB_paths}/eggnog_db_v5.0.2/eggnog.db
tar -xzvf ${DB_paths}/eggnog_db_v5.0.2/eggnog_proteins.dmnd.gz -C ${DB_paths}/eggnog_db_v5.0.2/eggnog_proteins.dmnd
