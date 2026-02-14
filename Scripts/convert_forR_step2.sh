#!/bin/bash

# Ensure that a parent argument is provided
if [ -z "$1" ]; then
  echo "No argument provided. Usage: $0 <parent>"
  exit 1
fi

parents=$1

# Ensure the output directory exists
output_dir="forR"
#mkdir -p "$output_dir"

awk -F "\t" -v OFS="\t" 'NR==1{print $1"_"$2,$9,$10}NR>1{if($7=="0/0") print $1"_"$2,$4,$3;else print $1"_"$2,$3,$4}' $output_dir/tmp.$parents.parents.MT.homo.vcf >$output_dir/$parents.parents.MT.forR
  
chrom=$(seq 1 25)
for chr in $chrom; do
  input_file="$output_dir/tmp.$parents.parents.$chr.homo.vcf"
  output_file="$output_dir/$parents.parents.$chr.forR"

  # Check if the input file exists
  if [ -f "$input_file" ]; then
    awk -F "\t" -v OFS="\t" 'NR==1{print $1"_"$2,$9,$10}NR>1{if($7=="0/0") print $1"_"$2,$4,$3;else print $1"_"$2,$3,$4}' "$input_file" >"$output_file"
  else
    echo "Warning: $input_file not found, skipping."
  fi
done  
echo "####### convert step1 finished ########"