#!/bin/bash

# Ensure that a parent argument is provided
if [ -z "$1" ]; then
  echo "No argument provided. Usage: $0 <parent>"
  exit 1
fi

parents=$1

# Ensure the output directory exists
output_dir="/forR"
mkdir -p "$output_dir"

zgrep -v "^##" homo/$parents.parents.MT.homo.vcf.gz|awk -F "\t" -v OFS="\t" 'NR==1{print $1,$2,$4,$5,$10,"TL_AD",$11,"TU_AD"}NR>1{split($10,TTTL,":");split($11,TTTU,":");print $1,$2,$4,$5,TTTL[1],TTTL[2],TTTU[1],TTTU[2]}'|awk -F "\t" -v OFS="\t" 'NR==1{print $0}NR>1{if(($5=="0/0" && $7=="1/1")||($7=="0/0" && $5=="1/1"))print $0}' |awk -F "\t" -v OFS="\t" 'NR==1{print $0,"TL_base","TU_base"}NR>1{if($5=="0/0")print $0,$3,$4;else print$0,$4,$3}' >$output_dir/tmp.$parents.parents.MT.homo.vcf
  
chrom=$(seq 1 25)
for chr in $chrom; do
  input_file="homo/$parents.parents.$chr.homo.vcf.gz"
  output_file="$output_dir/tmp.$parents.parents.$chr.homo.vcf"

  # Check if the input file exists
  if [ -f "$input_file" ]; then
    zgrep -v "^##" "$input_file"|awk -F "\t" -v OFS="\t" 'NR==1{print $1,$2,$4,$5,$10,"TL_AD",$11,"TU_AD"}NR>1{split($10,TTTL,":");split($11,TTTU,":");print $1,$2,$4,$5,TTTL[1],TTTL[2],TTTU[1],TTTU[2]}'|awk -F "\t" -v OFS="\t" 'NR==1{print $0}NR>1{if(($5=="0/0" && $7=="1/1")||($7=="0/0" && $5=="1/1"))print $0}' |awk -F "\t" -v OFS="\t" 'NR==1{print $0,"TL_base","TU_base"}NR>1{if($5=="0/0")print $0,$3,$4;else print$0,$4,$3}' >"$output_file"
  else
    echo "Warning: $input_file not found, skipping."
  fi
done  
echo "####### convert step1 finished ########"