#/bin/bash

# MT
./gatk-4.1.8.0 GenomicsDBImport -V TLCM.MT.vcf.gz -V TUAF1.MT.vcf.gz --genomicsdb-workspace-path UL.parents_chromMT  --intervals MT &


chrom=$(seq 1 25)
echo "start time: $(date)"
for chr in $chrom
do
        ./gatk-4.1.8.0 GenomicsDBImport -V TLCM.$chr.vcf.gz -V TUAF1.$chr.vcf.gz --genomicsdb-workspace-path UL.parents_chrom$chr --intervals $chr &
done
wait
echo "end time: $(date)"
echo "######## combine TUAF & TLCM GVCFs finished #######"
