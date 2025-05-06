library(ggplot2)
library(readr)
library(dplyr)

# Function to plot the distribution of -log10(p) values across chromosomes
args <- commandArgs(trailingOnly = TRUE)
input_file_path <- args[1]
output_dir <- args[2]
trait1 <- args[3]
trait2 <- args[4]

output_file_path_plot <- file.path(output_dir, paste0("LocalGeneticCorrelation_pval_", trait1, "_", trait2, ".png"))


df <- read.table(input_file_path, header = TRUE, stringsAsFactors = FALSE)
print(head(df))

df <- df %>% mutate(neg_log10_p = -log10(p))

chr_counts <- df %>%
  group_by(chr) %>%
  summarise(n = n()) %>%
  mutate(chr = as.factor(chr))

p <- ggplot(df, aes(x = factor(chr), y = neg_log10_p)) +
  geom_violin(fill = "skyblue", color = "black", scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, outlier.color = "red") +
  geom_text(data = chr_counts, aes(x = chr, y = max(df$neg_log10_p) + 0.5, label = n),
            inherit.aes = FALSE, size = 3, vjust = 0) +
  labs(title = "-log10(p) Distribution Across Chromosomes",
       subtitle = "Number of regions shown above violins",
       x = "Chromosome",
       y = "-log10(p)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave(output_file_path_plot, plot = p, width = 10, height = 6) 


## Make a table 

summary_table <- df %>%
  group_by(chr) %>%
  summarise(
    n_regions = n(),
    avg_p = mean(p, na.rm = TRUE),
    avg_h2_1 = mean(h2_1, na.rm = TRUE),
    avg_h2_2 = mean(h2_2, na.rm = TRUE)
  ) %>%
  arrange(as.numeric(chr))

# Save to a file
output_table_path <- file.path(output_dir, paste0("SummaryTable_", trait1, "_", trait2, ".tsv"))

if (file.exists(output_table_path)) {
  # If file exists, append without writing column names
  write.table(summary_table, output_table_path, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)
} else {
  # If file doesn't exist, create a new file with column names
  write.table(summary_table, output_table_path, sep = "\t", row.names = FALSE, col.names = TRUE, append = FALSE, quote = FALSE)
}





