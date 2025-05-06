library(data.table)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
trait <- args[1]
input_file <- args[2]
output_dir <- args[3]


p_value_cutoff <- c(1.00E-10, 5.00E-10, 1.00E-09, 5.00E-09, 1.00E-08, 5.00E-08, 1.00E-07, 5.00E-07, 1.00E-06, 5.00E-06, 1.00E-05, 5.00E-05, 1.00E-04, 5.00E-04, 1.00E-03, 5.00E-03, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1) 
k <- 4

for (ite in 1:k) {
    df <- fread(paste0(input_file, "/", trait, ".gwas.ite", ite, ".txt"))

    weights_df <- df %>% select(CHR, SNP, A1, A2)

    for (p_thr in p_value_cutoff) {
        # Create a column name that reflects the p-value threshold. You can format it as needed.
        colname <- paste0("P_", format(p_thr, scientific = TRUE))
        # Add the column: if the SNP's p-value is below the threshold, weight = beta; otherwise, weight = 0
        weights_df[, (colname) := ifelse(df$P <= p_thr, df$BETA, 0)]
    }
    outfile <- paste0(output_dir, "/", trait, ".P_and_T.ite", ite, ".txt")
    fwrite(weights_df, outfile, col.names=T, row.names=F, sep = "\t", quote=F)
}
