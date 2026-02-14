#/bin/bash

# MT
./gatk-4.1.8.0 GenotypeGVCFs -R Danio_rerio.GRCz11.fa -V gendb://UL.parents_chromMT -G StandardAnnotation -O UL.parents.MT.combine.vcf.gz &

chrom=$(seq 1 25)
echo "start time: $(date)"
for chr in $chrom
do
        ./gatk-4.1.8.0 GenotypeGVCFs -R Danio_rerio.GRCz11.fa -V gendb://UL.parents_chrom$chr -G StandardAnnotation -O UL.parents.$chr.combine.vcf.gz &
done
wait
echo "end time: $(date)"
echo "####### joint TUAF & TLCM genotyping finished ########"