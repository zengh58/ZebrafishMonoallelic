#!/bin/bash
SAMPLE=$1
echo "Modify cellsnp output: $(date)"
zgrep -v "^##" ${SAMPLE}/cellSNP.cells.vcf.gz |awk -F "\t" -v OFS="\t" 'NR==1{split($0,NNName,"\t")}NR>1{for(i=10;i<=NF;i=i+1){ if($i!=".:.:.:.:.:."){ split($i,DDDetail,":");{split($8,BBB,";")};print $1,$2,$4,$5,BBB[1],BBB[2],DDDetail[1],DDDetail[2],DDDetail[3],NNName[i] } }}' >${SAMPLE}.cell2snp
echo "Modify Done: $(date)"
#######################################
echo "Modify for R: $(date)"
awk -F "\t" -v OFS="\t" '{print $1,$2,$1"_"$2,$3,$4,$9-$8,$8,$7,$10}' ${SAMPLE}.cell2snp> ${SAMPLE}.cell2snp.forR
echo "Modify for R Done: $(date)"
##########################################
echo "Divide as Chromosome: $(date)"
mkdir "${SAMPLE}_cellsnp/"
awk -F "\t" -v OFS="\t" '{print > ("'"${SAMPLE}"'_cellsnp/chr." $1 ".cellsnp")}' "${SAMPLE}.cell2snp.forR"
echo "All done: $(date)"