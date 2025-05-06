trait=$1
PUMAS_exe_path=$2
k=$3
QC_input=$4
ld_path=$5
output=$6

Rscript $PUMAS_exe_path/PUMAS.subsampling.R \
--k $k \
--partitions 0.75,0.25 \
--trait_name $trait \
--gwas_path $QC_input \
--ld_path $ld_path \
--parallel \
--output_path $output/$trait/

