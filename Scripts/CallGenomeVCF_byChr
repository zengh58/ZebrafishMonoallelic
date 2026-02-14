#!/bin/bash
SM=$1
REF="Danio_rerio.GRCz11.fa"
outdir="4_gatk_byChr"
inputdir="3_markdup"

rm ${outdir}/${SM}/${SM}.log2
touch ${outdir}/${SM}/${SM}.log2

############################################################################
chroms1=$(seq 1 7)  
echo ${chroms1}
echo "1-7 start time: $(date)"
for chr in ${chroms1}
do
                echo ${chr}
                        ./gatk-4.1.8.0 HaplotypeCaller -R ${REF} -I ${inputdir}/${SM}.markdup.bam --intervals ${chr} --standard-min-confidence-threshold-for-calling 30 -A AlleleFraction -ERC GVCF -O ${outdir}/${SM}/${SM}.${chr}.vcf.gz &

                done 2>>${outdir}/${SM}/${SM}.log2
                wait
                echo "chrom 1-7 done"
                echo "1-7 end time: $(date)"
###########################################################################
chroms2=$(seq 8 16)
echo ${chroms2}
echo "8-16 start time: $(date)"
for chr in ${chroms2}
do
                echo ${chr}
                        ./gatk-4.1.8.0 HaplotypeCaller -R ${REF} -I ${inputdir}/${SM}.markdup.bam --intervals ${chr} --standard-min-confidence-threshold-for-calling 30 -A AlleleFraction -ERC GVCF -O ${outdir}/${SM}/${SM}.${chr}.vcf.gz &
                done 2>>${outdir}/${SM}/${SM}.log2
                wait
                echo "chrom 8-16 done"
                echo "8-16 end time: $(date)"
############################################################################
chroms3=$(seq 17 25)
echo ${chroms3}
echo "17-25 start time: $(date)"
for chr in ${chroms3}
do
                echo ${chr}
                        ./gatk-4.1.8.0 HaplotypeCaller -R ${REF} -I ${inputdir}/${SM}.markdup.bam --intervals ${chr} --standard-min-confidence-threshold-for-calling 30 -A AlleleFraction -ERC GVCF -O ${outdir}/${SM}/${SM}.${chr}.vcf.gz &
                done 2>>${outdir}/${SM}/${SM}.log2
                wait
                echo "17-25 end time: $(date)"
                echo "all chroms done"
