# Script: DESeq2-analyse + Volcano plot

library(DESeq2)
library(EnhancedVolcano)

counts <- read.csv("count_matrix_groot.csv", row.names = 1)
colnames(counts) <- paste0("ra", 1:8)
counts <- round(as.matrix(counts))

treatment <- c("Normal", "Normal", "Normal", "Normal", "RA", "RA", "RA", "RA")
treatment_table <- data.frame(treatment)
rownames(treatment_table) <- colnames(counts)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = treatment_table,
                              design = ~ treatment)
dds <- DESeq(dds)
resultaten <- results(dds)
write.table(resultaten, file = "Resultaten_RA_vs_Normal.csv", row.names = TRUE, col.names = TRUE)

png("VolcanoplotWC.png", width = 8, height = 10, units = "in", res = 500)
EnhancedVolcano(resultaten, lab = rownames(resultaten), x = 'log2FoldChange', y = 'padj')
dev.off()
