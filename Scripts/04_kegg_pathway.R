# KEGG (Scripts/04_kegg_pathway.R)
library(pathview)
gene_vector <- resultaten$log2FoldChange
names(gene_vector) <- rownames(resultaten)

pathview(gene.data = gene_vector, 
         pathway.id = "hsa05323", 
         species = "hsa", 
         gene.idtype = "SYMBOL", 
         limit = list(gene = 5))
