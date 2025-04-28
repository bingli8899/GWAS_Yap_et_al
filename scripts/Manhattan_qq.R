library(qqman)
library(dplyr)

snp_info <- read.table("data/sumstats/baldness.tsv", header = TRUE, sep = "\t")
sumstats <- read.table(gzfile("results/step1.sumstats.gz"), header = TRUE, sep = "\t", comment.char = "")

sumstats_chr <- left_join(sumstats, snp_info[, c("SNP", "CHR")], by = "SNP")

snp_info_clean <- snp_info[
  !is.na(snp_info$P_BOLT_LMM_INF) &
  is.finite(snp_info$P_BOLT_LMM_INF) &
  snp_info$P_BOLT_LMM_INF > 0 &
  snp_info$P_BOLT_LMM_INF <= 1, 
]

png("results/plots/manhattan_plot.png", width = 1200, height = 600)
manhattan(snp_info, chr = "CHR", bp = "BP", snp = "SNP", p = "P_BOLT_LMM_INF",
          genomewideline = -log10(5e-8), suggestiveline = -log10(1e-5),
          main = "Manhattan Plot")
dev.off()

png("results/plots/qq_plot.png", width = 600, height = 600)
qq(snp_info$P_BOLT_LMM_INF, main = "QQ Plot")
dev.off()


