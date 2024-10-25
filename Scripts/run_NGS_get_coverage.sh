#!/usr/bin/env bash

set -euo pipefail

# Data directory


if [[ -z "$1" ]]; then
  log "Must supply the data path.\n"
  exit 1
else
  LOCAL_DATA_DIR="$1"
fi

# Check if the directory exists
if [ -d "$LOCAL_DATA_DIR" ]; then
    echo "Directory exits"
    
else
    echo "Error: $LOCAL_DATA_DIR is not a valid directory."
    exit 1
fi

# the target file
input_suffix=".cram"


for file1 in "${LOCAL_DATA_DIR}/Data/"*"${input_suffix}"; do
  echo "Process the data: ${file1}"

  sample="$(echo "$file1" | sed 's/\.cram//' | sed 's/Data/Results/')" 
  echo $sample
  
  echo "Step 1: using Bedtools to get the coverage data: # of reads overlapping the position"
  bedtools genomecov -ibam ${file1} -bg > ${sample}.bedtools_cov.txt 

  echo "Step 2: using Bedtools to get the coverage histogram file: percent of reads overlapping the position"
  bedtools genomecov -ibam ${file1} > ${sample}.bedtools_hist.txt 

  echo "Step 3a:  Extract the overall coverage stats from the coverage file"
  awk '{sum+=$4; count++} END {print "Total number of reads across all read-covered positions: ", sum, "\nTotal number of read-covered positions: ", count}' ${sample}.bedtools_cov.txt > ${sample}.bedtools_cov_sum_stat.txt

  echo "Step 3b: Total covered region: "
  awk '{ if ($4 > 0) { sum += $3 - $2 } } END { print "Total Length of Read Covered Genome Regions (bp):", sum,"\nThe proportion of GRCh38 has coverage (%):", sum*100/(3.2*10^9)}' ${sample}.bedtools_cov.txt >> ${sample}.bedtools_cov_sum_stat.txt

done