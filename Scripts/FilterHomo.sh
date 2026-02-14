#!/bin/bash

# Ensure that a parent argument is provided
if [ -z "$1" ]; then
  echo "No argument provided. Usage: $0 <parent>"
  exit 1
fi

parents=$1

# Ensure the output directory exists
output_dir="homo"
mkdir -p "$output_dir"

zcat modified/$parents.parents.MT.mod.vcf.gz| awk -F "\t" -v OFS="\t" 'NR<=1033{print $0} NR>1033{split($10,TTTL,":");split($11,TTTU,":"); if((TTTL[1]=="0/0" && TTTU[1]=="1/1")||(TTTU[1]=="0/0" && TTTL[1]=="1/1")) print $0}' |gzip >$output_dir/$parents.parents.MT.homo.vcf.gz
  
chrom=$(seq 1 25)
for chr in $chrom; do
  input_file="modified/$parents.parents.$chr.mod.vcf.gz"
  output_file="$output_dir/$parents.parents.$chr.homo.vcf.gz"

  # Check if the input file exists
  if [ -f "$input_file" ]; then
    zcat "$input_file"|awk -F "\t" -v OFS="\t" 'NR<=1033{print $0}NR>1033{split($10,TTTL,":");split($11,TTTU,":");if((TTTL[1]=="0/0" && TTTU[1]=="1/1")||(TTTU[1]=="0/0" && TTTL[1]=="1/1"))print $0}' |gzip > "$output_file"
  else
    echo "Warning: $input_file not found, skipping."
  fi
done  
echo "####### filter homosite finished ########"