# Script: Sorteren, indexeren en gen-tellingen

library(Rsamtools)
samples <- paste0("ra", 1:8)
lapply(samples, function(s) sortBam(file = paste0(s, ".BAM"), destination = paste0(s, ".sorted")))

library(readr)
library(dplyr)
library(Rsubread)

gff <- read_tsv("GCF_000001405.25_GRCh37.p13_genomic.gtf.gz", comment = "#", col_names = FALSE)
colnames(gff) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
gff_gene <- gff %>% filter(type == "gene")
gff_gene$type <- "exon"

allsamples <- paste0("ra", 1:8, ".bam")

count_matrix <- featureCounts(files = allsamples,
                              annot.ext = "GCF_000001405.25_GRCh37.p13_genomic.gtf.gz",
                              isPairedEnd = TRUE,
                              isGTFAnnotationFile = TRUE,
                              GTF.attrType = "gene_id",
                              useMetaFeatures = TRUE)

write.csv(count_matrix$counts, "count_matrix_groot.csv")
