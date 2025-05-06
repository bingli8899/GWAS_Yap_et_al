library(data.table)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
trait <- args[1]
input_file <- args[2]
pruned_file <- args[3]
output_dir <- args[4]

gwas.tmp <- fread(paste0(input_file, trait,".txt.gz"))

prune <- fread(paste0(pruned_file, trait, ".prune.in"), header=F)$V1

gwas <- gwas.tmp %>% filter(SNP %in% prune)

nrow(gwas)

fwrite(as.data.frame(gwas), paste0(output_dir, trait, ".txt.gz"), col.names=T, row.names=F, sep="\t", quote=F)
