# ZebrafishMonoallelic

In this study, we present a comprehensive single-cell atlas of allelic expression in zebrafish and systematically characterize allele-specific expression patterns using a maximum likelihood estimation (MLE) framework.

## PROCEDURE
### 1. **Parental Whole Genome Sequencing (WGS) variant calling**

* Read alignment
```bash
trim_galore --phred33 -q 25 --length 35 --stringency 3 --fastqc --paired --max_n 3 -o 1_trim_galore TUAF1_combined_R1.fastq.gz TUAF1_combined_R2.fastq.gz

bwa mem -t 8 -M genome/Danio_rerio.GRCz11.fa 1_trim_galore/TUAF1_combined_R1_val_1.fq.gz 1_trim_galore/TUAF1_combined_R2_val_2.fq.gz -R "@RG\tID:TUAF1\tLB:TUAF1\tSM:TUAF1\tPL:ILLUMINA" >2_bwa/TUAF1.sam

samtools view -b -u 2_bwa/TUAF1.sam|samtools sort -@ 8 -o 2_bwa/TUAF1.sorted.bam 

samtools index -@ 15 2_bwa/TUAF1.sorted.bam
```
* Variant calling using GATK
```bash
./gatk-4.1.8.0 MarkDuplicates -I 2_bwa/TUAF1.sorted.bam -O 3_markdup/TUAF1.markdup.bam -M 3_markdup/TUAF1.metrics --CREATE_INDEX &

# chr 1-25
./callGenomeVCF_byChr.sh TUAF1
# chr MT
./gatk-4.1.8.0 HaplotypeCaller -R Danio_rerio.GRCz11.fa -I 3_markdup/TUAF1.markdup.bam --intervals MT --standard-min-confidence-threshold-for-calling 30 -A AlleleFraction -ERC GVCF -O 4_gatk_byChr/TUAF1/TUAF1.MT.vcf.gz

./combine_UL_vcf_step1.sh
./combine_UL_vcf_step2.sh
```
* Retention of heterozygous sites
```bash
# replace | with /
./replace.sh UL

# Filter
./FilterHomo.sh UL

# Merge multiple chromosome vcf
./Merge.sh

# Filter GQ & SNP
./gatk-4.1.8.0 VariantFiltration \
-V homo/UL.parents.allhomo.vcf.gz \
-filter "GQ < 20" --filter-name "lowGQ" \
-O UL.parents.filter.homo.vcf.gz
```
* Conversion of variant data into R-compatible format
```bash
./convert_forR_step1.sh UL
./convert_forR_step2.sh UL

cat forR/UL.parents.*.forR > forR/UL.parents.all.forR
```

### 2. **Single-cell RNA-seq data processing**
All scRNA-seq data were processed following the standard pipeline provided by [10x genenomics](https://www.10xgenomics.com/support/cn/software/cell-ranger/latest).

```bash
 cellranger count --id=LU1_cellranger \
--fastqs=Sample_LU1 \
--sample=LU1 \
--transcriptome=ref/Danio.rerio_genome \
--localcores=15 \
--localmem=64 \
--include-introns=true \
--create-bam=true
```

### 3. **Cell-level SNV expression profiling with cellSNP-lite tools**
1. Cellsnp-lite is a C/C++ tool for efficient genotyping bi-allelic SNPs on single cells. Please refer to the official [cellsnp-lite](https://cellsnp-lite.readthedocs.io/en/latest/) documentation for detailed usage.

   ```bash
   cellsnp-lite -s LU1.possorted_genome_bam.bam \
   -b LU1.barcodes.tsv.gz \
   -O LU1 -R UL.parents.filter.homo.vcf.gz \
   -p 120 --minMAF 0.1 --minCOUNT 20 \
   --gzip --genotype &
   ```

2. Conversion of variant data into R-compatible format.
   
   ```bash
   ./Modifycellsnp.sh LU1
   ```
   
### 4. **Maximum Likelihood Estimation Framework**
* MLE_code.R


