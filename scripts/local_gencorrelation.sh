#!/bin/bash

sumstats1=$1
sumstats2=$2
N1=$3
N2=$4
bfile_path=$5
partition_path=$6
output_path=$7
python_path=$8
supergnova_path=$9

"$python_path" "$supergnova_path/supergnova.py" \
    "$sumstats1" "$sumstats2" \
    --N1 "$N1" \
    --N2 "$N2" \
    --bfile "$bfile_path" \
    --partition "$partition_path" \
    --out "$output_path"
