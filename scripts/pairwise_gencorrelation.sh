#!/bin/bash

# Arguments
ldsc_path=$1          # Path to ldsc-2.0.1
sumstats1=$2          # Full path to first sumstats.gz file
sumstats2=$3          # Full path to second sumstats.gz file
ref_path=$4           # Path to reference folder
output_path=$5        # Output file prefix (no extension)
python_path=$6        # Full path to python executable
trait1=$7         # Name of the first trait --> Should always be baldness in our case
trait2=$8         # Name of the second trait --> Should always be baldness in our case

# Run ldsc.py for genetic correlation
"$python_path" "$ldsc_path/ldsc.py" \
  --rg "$sumstats1","$sumstats2" \
  --ref-ld-chr "$ref_path/eur_w_ld_chr/" \
  --w-ld-chr "$ref_path/eur_w_ld_chr/" \
  --out "$output_path/${trait1}_${trait2}.step4"