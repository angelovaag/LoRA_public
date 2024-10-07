#!/bin/bash
#source: https://github.com/EnvGen/metagenomics-workshop/blob/master/in-house/prokkagff2gtf.sh
#Usage: prokkagff2gtf.sh <prodigal gff file>

infile=$1

if [ "$infile" == "" ] ; then
    echo "Usage: prodigalgff2gtf.sh <prodigal gff file>"
    exit 0
fi

#old version
#grep -v "#" $infile | grep "ID=" | cut -f1 -d ';' | sed 's/ID=//g' | cut -f1,4,5,7,9 | \
# awk -v OFS='\t' '{print $1,"Prodigal","CDS",$2,$3,".",$4,".","gene_id " $5}'

#my version
grep -v "#" $infile | grep "ID=" | cut -f1 -d ';' | sed 's/ID=//g' | cut -f1,4,5,7,9 |  \
    awk -v OFS='\t' '{print $1,"Prodigal","CDS",$2,$3,".",$4,".","gene_id " $1 ":" $5}' |sed 's/:.*_/_/g'