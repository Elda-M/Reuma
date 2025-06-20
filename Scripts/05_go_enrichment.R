# GO analyse (Scripts/05_go_enrichment.R)
library(org.Hs.eg.db)
library(GO.db)
library(goseq)
library(clusterProfiler)

converted <- bitr(rownames(resultaten), fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = org.Hs.eg.db)
resultaten <- as.data.frame(resultaten)
resultaten$symbol <- rownames(resultaten)
resultaten <- dplyr::left_join(resultaten, converted, by = c("symbol" = "SYMBOL"))

all_genes <- as.integer(resultaten$padj < 0.05 & abs(resultaten$log2FoldChange) > 1)
names(all_genes) <- resultaten$ENSEMBL

library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_lengths <- getBM(attributes = c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = names(all_genes), mart = ensembl)

gene_lengths_unique <- gene_lengths %>% group_by(ensembl_gene_id) %>% summarise(length = round(mean(transcript_length, na.rm = TRUE)))
bias.data <- gene_lengths_unique$length
names(bias.data) <- gene_lengths_unique$ensembl_gene_id

gemeenschappelijke_genen <- intersect(names(all_genes), names(bias.data))
all_genes_clean <- all_genes[gemeenschappelijke_genen]
bias.data_clean <- bias.data[gemeenschappelijke_genen]

pwf <- nullp(DEgenes = all_genes_clean, 
             genome = "hg38", 
             id = "ENSEMBL", 
             bias.data = bias.data_clean)

png("pwf_plot.png", width = 800, height = 600)
plotPWF(pwf)
dev.off()

gene2go <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(pwf), keytype = "ENSEMBL", columns = c("GO", "ONTOLOGY"))
gene2go <- gene2go[gene2go$ONTOLOGY == "BP", ]
gene2cat <- split(gene2go$GO, gene2go$ENSEMBL)

GO_results <- goseq(pwf, gene2cat = gene2cat, method = "Hypergeometric", use_genes_without_cat = TRUE)

library(ggplot2)
png("GO_resultaten_plot.png", width = 1000, height = 600)
GO_results %>% top_n(10, wt = -over_represented_pvalue) %>% mutate(hitsPerc = numDEInCat * 100 / numInCat) %>% ggplot(aes(x = hitsPerc, y = term, colour = over_represented_pvalue, size = numDEInCat)) + geom_point() + expand_limits(x = 0) + labs(x = "Hits (%)", y = "GO term", colour = "p value", size = "Count")
dev.off()

write.csv(GO_results, file = "GO_resultaten.csv", row.names = FALSE)
