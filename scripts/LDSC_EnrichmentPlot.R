source("scripts/CellType.R") 
library(ggplot2)     
library(ggthemes) 
library(flextable)



results <- read.table("results/step3.results", header=TRUE)
cell_rows <- grep("L2_1$", results$Category)
celltype_table <- results[cell_rows, c("Category", "Enrichment", "Enrichment_std_error", "Enrichment_p")]
colnames(celltype_table) <- c("Cell_Type", "Enrichment", "SE", "P_value")
celltype_table <- celltype_table[order(-celltype_table$Enrichment), ]
celltype_table$Enrichment <- round(celltype_table$Enrichment, 2)
celltype_table$SE <- round(celltype_table$SE, 3)
celltype_table$P_value <- signif(celltype_table$P_value, 3)
print(celltype_table)
ft <- flextable(celltype_table)
save_as_docx(ft, path = "celltype_enrichment_table.docx")


## input:
# Directory to the LDSC heritability partitioning results (change this):
InDir="results"
# LDSC partition heritability results (change this):
h2.result="step3.results"

## output:

# Directory to store the enrichment plot (change this):
OutDir="results/plots"
# enrichment plot file name (change this):
out.fn="enrichment_plot"

# no need to change below unless you know what you are doing:
Colors.HF66 = c(rep("#FB8575",18),"#D17493","#8D719B",rep("#446B85",8),rep("#229E9C",2),"#5DCD8F",rep("#CAF270",3),rep("#DAC14D",5),
                rep("#D79347",12),rep("#C06B4D",2),"#677C36",rep("#177374",3),rep("#6D5271",3),rep("#6389AE",2),rep("#34C5CD",3),"#62FCC0")


# tissue in the Cells.HF so use the following chunk of codes:
Data = read.table(paste0(InDir, "/", h2.result), header=T, stringsAsFactors=F)
Data = data.frame(Tissue=factor(Cells.HF[-c(32,62:64)], levels=unique(Cells.HF[-c(32,62:64)])), Pvalue=as.numeric(Data[c(55:120),7]))

# can adjust the size of the output png file:
png(paste0(OutDir, "/", out.fn, ".png"), width=800, height=1200, type="cairo")
gg <- ggplot(Data, aes(x=Tissue, y=-log10(Pvalue)))
gg <- gg + geom_bar(stat="identity", fill=rev(Colors.HF66)) + coord_flip() + ggtitle("All")
gg <- gg + scale_y_continuous(expand=c(0,0))
gg <- gg + theme_minimal(base_family="Helvetica")
gg <- gg + theme(panel.grid.major.y=element_blank(), panel.grid.minor=element_blank())
gg <- gg + theme(plot.title=element_text(size=22, face="bold"), axis.ticks=element_blank(), axis.text=element_text(size=18), axis.title=element_text(size=22))
gg <- gg + labs(x=NULL, y=expression("-log"[10]*"(p-value)"))
gg <- gg + geom_hline(yintercept = -log10(0.05/(66)), linetype=5, alpha=0.5)
print(gg)
dev.off()

q(save = "no")

