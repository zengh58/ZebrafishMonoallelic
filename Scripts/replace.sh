#!/bin/bash

# Ensure that a parent argument is provided
if [ -z "$1" ]; then
  echo "No argument provided. Usage: $0 <parent>"
  exit 1
fi

parents=$1
chrom=$(seq 1 25)

# Ensure the output directory exists
output_dir="modified"
mkdir -p "$output_dir"

for chr in $chrom; do
  input_file="$parents.parents.$chr.combine.vcf.gz"
  output_file="$output_dir/$parents.parents.$chr.mod.vcf.gz"

  # Check if the input file exists
  if [ -f "$input_file" ]; then
    zcat "$input_file" | awk 'BEGIN {FS=OFS="\t"} {sub(/\|/, "/", $10); sub(/\|/, "/", $11); print}' | gzip > "$output_file"
  else
    echo "Warning: $input_file not found, skipping."
  fi
done  
echo "####### replace finished ########"