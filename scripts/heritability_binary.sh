#!/bin/bash

# Arguments
input_path=$1        # Path to folder containing the .txt.gz file
ldsc_path=$2         # Path to ldsc-2.0.1 directory
Ncas=$3              # Number of cases
Ncon=$4              # Number of controls
ref_path=$5          # Path to reference folder
output_dir=$6        # Output directory
SNP=$7               # SNP column name
ALLELE1=$8           # Allele1 column name
ALLELE0=$9           # Allele0 column name
P=${10}              # P-value column name
sumstats=${11}       # Signed sumstats column (e.g., LOG_OR)
data_name=${12}      # Data file name (without extension)
python_path=${13}    # Full path to python executable

# Make sure output directory exists
mkdir -p "$output_dir"

# Run munge_sumstats.py
"$python_path" "$ldsc_path/munge_sumstats.py" \
  --sumstats "$input_path/$data_name.tsv" \
  --N-cas "$Ncas" \
  --N-con "$Ncon" \
  --signed-sumstats "$sumstats",0 \
  --snp "$SNP" \
  --a1 "$ALLELE1" \
  --a2 "$ALLELE0" \
  --p "$P" \
  --out "$output_dir/$data_name.step1" \
  --merge-alleles "$ref_path/w_hm3.snplist"
