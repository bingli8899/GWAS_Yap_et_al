#!/bin/bash

# bash script for heritability analysis using ldsc

input_path=$1
ldsc_path=$2
N=$3
ref_path=$4
output_dir=$5

SNP=$6
ALLELE1=$7
ALLELE0=$8
P=$9
sumstats=${10}
data_name=${11}
python_path=${12}

# conda activate ldsc

mkdir -p "$output_dir"

"$python_path" "$ldsc_path/munge_sumstats.py" \
    --sumstats "$input_path/$data_name.tsv" \
    --N "$N" \
    --signed-sumstats "$sumstats,0" \
    --snp "$SNP" \
    --a1 "$ALLELE1" \
    --a2 "$ALLELE0" \
    --p "$P" \
    --merge-alleles "$ref_path/w_hm3.snplist" \
    --out "$output_dir/$data_name.step1"

if [ "$data_name" == "baldness" ]; then
  echo "Running on baldness dataset..."

  "$python_path" "$ldsc_path/ldsc.py" \
    --h2 "$output_dir/$data_name.step1.sumstats.gz" \
    --ref-ld-chr "$ref_path/eur_w_ld_chr/" \
    --w-ld-chr "$ref_path/eur_w_ld_chr/" \
    --out "$output_dir/$data_name.step2"

  "$python_path" "$ldsc_path/ldsc.py" \
    --h2 "$output_dir/$data_name.step1.sumstats.gz" \
    --ref-ld-chr "$ref_path/Baseline/baseline.","$ref_path/GenoSkylinePlus/GSplus_Tier3_1KGphase3." \
    --w-ld-chr "$ref_path/weights/weights.hm3_noMHC." \
    --frqfile-chr "$ref_path/genotype/1000G.EUR.QC." \
    --overlap-annot \
    --out "$output_dir/$data_name.step3"
fi

# conda deactivate
