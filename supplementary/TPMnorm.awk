#!/bin/bash

## TMP normalization for preGENES


gtfFILE=${1}


cut -f4,5,9 ${gtfFILE} | sed 's/gene_id //g' | gawk '{print $3,$2-$1+1}' | tr ' ' '\t' | sed '1s/^/gene\tlength\n/' > {output.tmp_lengths} 
join {output.tmp_lengths} {input.verse_cds} -t $'\t' > {output.tmp_counts}
grep \"_\" {output.tmp_counts}   | awk -v OFS=\"\t\" '{$4 = sprintf(\"%0.0f\", $3*150/$2)}1' | sort | sed -e '1s/^/#name\tlen\treads\tcov\n/' > {output.tmp_coverage}                               ##--> do we even need this? <--##
grep \"_\" {output.tmp_coverage} | awk -v OFS=\"\t\" '{$5 = sprintf(\"%0.0f\", $3*1e3/$2)}1' | sort | sed -e '1s/^/#name\tlen\treads\tcov\tRPK\n/' > {output.tmp_rpk}                         ##--> updatad code <--
  awk -v OFS=\"\t\" 'NR==FNR{sum+= $5; next} FNR==1{print $0,"iTPM"; next} {printf("%s %0.0f\n", $0,$5*1e6/sum)}'  {output.tmp_rpk} {output.tmp_rpk} |sed -e 's/\\s/\\t/g' > {output.ABUNtab}     ##--> calc iTPM now <--
# awk -v OFS=\"\t\" 'NR==FNR{sum+= $5; next} FNR==1{print $0,"iTPM"; next} {printf("%s %0.0f\n", $0,$5*1e6/sum)}}'  {output.tmp_rpk} {output.tmp_rpk} |sed -e 's/\\s/\\t/g' > {output.ABUNtab}     ##--> calc iTPM now <--