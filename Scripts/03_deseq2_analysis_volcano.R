# DESeq2 analyse + volcano (Scripts/03_deseq2_analysis_volcano.R)
count_matrix_groot <- read.table("count_matrix.txt")
write.csv(count_matrix_groot, "count_matrix_groot.csv")

counts <- read.csv("count_matrix_groot.csv", row.names = 1)
colnames(counts) <- paste0("ra", 1:8)
counts <- round(as.matrix(counts))

treatment <- c("Normal", "Normal", "Normal", "Normal", "RA", "RA", "RA", "RA")
treatment_table <- data.frame(treatment)
rownames(treatment_table) <- colnames(counts)

library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = treatment_table, design = ~ treatment)
dds <- DESeq(dds)
resultaten <- results(dds)

write.table(resultaten, file = "Resultaten_RA_vs_Normal.csv", row.names = TRUE, col.names = TRUE)

library(EnhancedVolcano)
png("VolcanoplotWC.png", width = 8, height = 10, units = "in", res = 500)
EnhancedVolcano(resultaten, lab = rownames(resultaten), x = 'log2FoldChange', y = 'padj')
dev.off()