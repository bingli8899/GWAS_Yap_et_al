#!/bin/bash
out_dir=$1
bfile=$2
snp_lst=$3
plink_folder=$4
echo "Running pruning with output prefix $out_dir"

# Perform LD pruning using a window of 100 variants,
# shifting the window by 5 variants, and an r2 threshold of 0.1.
# Note: This command will output two files:
#   $2.prune.in  -> SNPs that passed the pruning (i.e., the independent set)
#   $2.prune.out -> SNPs that were removed due to high LD

$plink_folder/plink2 \
  --bfile $bfile \
  --extract $snp_lst \
  --indep-pairwise 1000 5 0.1 \
  --out $out_dir
  