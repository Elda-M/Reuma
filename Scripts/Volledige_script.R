setwd("C:/Users/Elda Mehmedagic/School")
list.files("C:/Users/Elda Mehmedagic/School")

install.packages("R.utils")  # Alleen de eerste keer nodig
library(R.utils)

# Vervang de bestandsnaam hieronder met je eigen zip-bestand
unzip("Data_RA_raw.zip", exdir = "Data_RA_raw") #Hiermee worden de bestanden uitgepakt in een submap 'ethanol_data'

install.packages('BiocManager')

BiocManager::install('Rsubread')
library(Rsubread)

buildindex(
  basename = 'ref_human',
  reference = 'GCF_000001405.40_GRCh38.p14_genomic.fna',
  memory = 7000,
  indexSplit = TRUE)

# Mapping RA-samples handmatig (paired-end FASTQ)

# RA monsters
align.ra1 <- align(index = "ref_human",
                   readfile1 = "SRR4785819_1_subset40k.fastq",
                   readfile2 = "SRR4785819_2_subset40k.fastq",
                   output_file = "ra1.BAM")

align.ra2 <- align(index = "ref_human",
                   readfile1 = "SRR4785820_1_subset40k.fastq",
                   readfile2 = "SRR4785820_2_subset40k.fastq",
                   output_file = "ra2.BAM")

align.ra3 <- align(index = "ref_human",
                   readfile1 = "SRR4785828_1_subset40k.fastq",
                   readfile2 = "SRR4785828_2_subset40k.fastq",
                   output_file = "ra3.BAM")

align.ra4 <- align(index = "ref_human",
                   readfile1 = "SRR4785831_1_subset40k.fastq",
                   readfile2 = "SRR4785831_2_subset40k.fastq",
                   output_file = "ra4.BAM")

align.ra5 <- align(index = "ref_human",
                   readfile1 = "SRR4785979_1_subset40k.fastq",
                   readfile2 = "SRR4785979_2_subset40k.fastq",
                   output_file = "ra5.BAM")

align.ra6 <- align(index = "ref_human",
                   readfile1 = "SRR4785980_1_subset40k.fastq",
                   readfile2 = "SRR4785980_2_subset40k.fastq",
                   output_file = "ra6.BAM")

align.ra7 <- align(index = "ref_human",
                   readfile1 = "SRR4785986_1_subset40k.fastq",
                   readfile2 = "SRR4785986_2_subset40k.fastq",
                   output_file = "ra7.BAM")

align.ra8 <- align(index = "ref_human",
                   readfile1 = "SRR4785988_1_subset40k.fastq",
                   readfile2 = "SRR4785988_2_subset40k.fastq",
                   output_file = "ra8.BAM")

# Laad Rsamtools voor sorteren en indexeren
library(Rsamtools)

# Bestandsnamen van de monsters
samples <- paste0("ra", 1:8)

# Voor elk monster: sorteer en indexeer de BAM-file
# Sorteer BAM-bestanden
lapply(samples, function(s) {
  sortBam(file = paste0(s, ".BAM"), destination = paste0(s, ".sorted"))
})

library(readr)
library(dplyr)
library(Rsamtools)
library(Rsubread)

# Inlezen en filteren van GFF3-bestand
gff <- read_tsv("GCF_000001405.25_GRCh37.p13_genomic.gtf.gz", comment = "#", col_names = FALSE)

# Kolomnamen toevoegen
colnames(gff) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

# Alleen genregels selecteren
gff_gene <- gff %>% filter(type == "gene")

# 'type' aanpassen naar 'exon' zodat featureCounts het accepteert
gff_gene$type <- "exon"

# Je definieert een vector met namen van BAM-bestanden. Elke BAM bevat reads van een RNA-seq-experiment (bijv. behandeld vs. controle).

allsamples <- paste0("ra", 1:8, ".bam")

count_matrix <- featureCounts(
  files = allsamples,
  annot.ext = "GCF_000001405.25_GRCh37.p13_genomic.gtf.gz",
  isPairedEnd = TRUE,
  isGTFAnnotationFile = TRUE,
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE)

head(count_matrix$annotation)
head(count_matrix$counts)

View(count_matrix$counts)

 
# Bekijk eerst de structuur van het object

count_matrix_groot <- read.table("count_matrix.txt")

#Nu de matrix klaar is, sla je deze op. Zo kun je hem in Werkcollege 3 direct opnieuw inlezen.
write.csv(count_matrix_groot, "count_matrix_groot.csv")

#Bekijk de eerste paar rijen.
head(count_matrix_groot)

# WERKCOLLEGE 3 Statistiek en analyse

# Lees het bestand uit werkcollege 2 in
counts <- read.csv("count_matrix_groot.csv", row.names = 1)

# Behandelingstabel maken

colnames(counts) <- paste0("ra", 1:8)  # of juiste sample-namen
counts <- round(as.matrix(counts))    # 🔧 dit is de fix

treatment <- c("Normal", "Normal", "Normal", "Normal", "RA", "RA", "RA", "RA")
treatment_table <- data.frame(treatment)
rownames(treatment_table) <- colnames(counts)


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2", "KEGGREST"))


library(DESeq2)
library(KEGGREST)

# Statistiek

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = treatment_table,
                              design = ~ treatment)

dds <- DESeq(dds)
resultaten <- results(dds)

# Bekijk de top van de resultaten
head(resultaten)

# Optioneel: resultaten opslaan
write.table(resultaten, file = "Resultaten_RA_vs_Normal.csv", row.names = TRUE, col.names = TRUE)

# Samenvatting significante genen
# stap 1: Hoeveel genen zijn er écht veranderd? Hier tellen we hoeveel genen er significant op- of neer-gereguleerd zijn.

sum(resultaten$padj < 0.05 & resultaten$log2FoldChange > 1, na.rm = TRUE)
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange < -1, na.rm = TRUE)

# Stap 2: Welke genen springen eruit? Nu sorteren we het resultaat om te kijken naar de opvallendste genen
hoogste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = TRUE), ]
laagste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = FALSE), ]
laagste_p_waarde <- resultaten[order(resultaten$padj, decreasing = FALSE), ]

# Bekijk nu welke genen het belangrijkst zijn volgens de analyse.
head(laagste_p_waarde)

#Visualisatie
# Een volcano plot laat zien welke genen significant veranderen in expressie. Op de x-as staat de log2 fold change en op de y-as de -log10 van de aangepaste p-waarde (padj). Alleen genen met padj < 0,05 en log2fc van > 1 of < -1 worden gelabeled.
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  BiocManager::install("EnhancedVolcano")
}
library(EnhancedVolcano)

EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'padj')

#  Alternatief (alle genen zichtbaar, geen p-waarde cutoff):
EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0)

# Sla het figuur op.
dev.copy(png, 'VolcanoplotWC.png', 
         width = 8,
         height = 10,
         units = 'in',
         res = 500)
dev.off()

# KEGG Pathway-analyse

if (!requireNamespace("pathview", quietly = TRUE)) {
  BiocManager::install("pathview")
}
library(pathview)

# Visualiseer een KEGG-pathway
library(pathview)

# Stel dat je genenlijst een named vector is zoals:
gene_vector <- resultaten$log2FoldChange
names(gene_vector) <- rownames(resultaten)

# Visualiseer het KEGG-pathway voor RA
pathview(
  gene.data  = gene_vector,
  pathway.id = "hsa05323",       # Rheumatoid arthritis pathway
  species    = "hsa",            # Human
  gene.idtype = "SYMBOL",       # of "ENTREZ", afhankelijk van je gen-ID's
  limit = list(gene = 5)
)


# HET VERVOLG
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("goseq", "org.Hs.eg.db", "GO.db"))
BiocManager::install("clusterProfiler")
library(GO.db)
library(org.Hs.eg.db)
library(goseq)
library(clusterProfiler)

# Converteer SYMBOL → ENSEMBL
converted <- bitr(rownames(resultaten),
                  fromType = "SYMBOL",
                  toType   = "ENSEMBL",
                  OrgDb    = org.Hs.eg.db)

# Zet symbol als kolom en join de data

resultaten <- as.data.frame(resultaten)

resultaten$symbol <- rownames(resultaten)

resultaten <- dplyr::left_join(
  resultaten,
  converted,
  by = c("symbol" = "SYMBOL")
)

# Zorg dat resultaten$padj en resultaten$log2FoldChange bestaan: Maak een binaire vector van significante genen
all_genes <- as.integer(resultaten$padj < 0.05 & abs(resultaten$log2FoldChange) > 1)
names(all_genes) <- resultaten$ENSEMBL

head(names(all_genes))     # geen NULL meer

# bereken/match de genlengtes

if (!requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}

# Genlengtes ophalen via biomaRt
library(biomaRt)

# Kies het Ensembl-biomart en dataset voor mens
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Haal de transcriptlengtes op voor de genen in jouw lijst (names(all_genes))
gene_lengths <- getBM(
  attributes = c("ensembl_gene_id", "transcript_length"),
  filters    = "ensembl_gene_id",
  values     = names(all_genes),
  mart       = ensembl
)


# maak een 	een object dat verwijst naar de juiste Ensembl-database en species
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


# Zet de NCBI IDs van gene2cat_table om
ncbi_ids <- names(gene2cat_table)

mapping <- getBM(
  attributes = c("entrezgene_id", "ensembl_gene_id"),
  filters    = "entrezgene_id",
  values     = ncbi_ids,
  mart       = mart
)


# Genlengtes samenvatten en omzetten

library(dplyr)

head(gene_)

lengths_avg <- gene_lengths %>%
  group_by(ensembl_gene_id) %>%
  summarise(length = round(mean(transcript_length, na.rm = TRUE)))

bias.data <- lengths_avg$length
names(bias.data) <- lengths_avg$ensembl_gene_id

# Filter alleen genen die zowel in all_genes als in bias.data zitten
gemeenschappelijke_genen <- intersect(names(all_genes), names(bias.data))

# Subset beide objecten
all_genes_filtered <- all_genes[gemeenschappelijke_genen]
bias.data_filtered <- bias.data[gemeenschappelijke_genen]

# Bepaal alleen genen met géén NA in beide vectoren
gemeenschappelijke_genen <- names(all_genes_filtered)[!is.na(all_genes_filtered)]

# Subset beide vectoren
all_genes_clean     <- all_genes_filtered[gemeenschappelijke_genen]
bias.data_clean     <- bias.data_filtered[gemeenschappelijke_genen]

sum(is.na(all_genes_clean))           # 0
sum(is.na(bias.data_clean))           # 0
identical(names(all_genes_clean), names(bias.data_clean))  # TRUE


pwf <- nullp(
  DEgenes   = all_genes_clean,
  genome    = "hg38",
  id        = "ENSEMBL",
  bias.data = bias.data_clean
)


# Sla direct de figuur op in een bestand
png("pwf_plot.png", width = 800, height = 600)
plotPWF(pwf)
dev.off()

# GO-analyse 
# Nodige libraries
library(org.Hs.eg.db)
library(GO.db)

# 1. Gen-ID's (ENSEMBL)
ensembl_ids <- rownames(pwf)

# 2. Haal GO mappings op
gene2go <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = ensembl_ids,
  keytype = "ENSEMBL",
  columns = c("GO", "ONTOLOGY")
)

# 3. Filter alleen GO:BP
gene2go <- gene2go[gene2go$ONTOLOGY == "BP", ]

# 4. Maak gene2cat object
gene2cat <- split(gene2go$GO, gene2go$ENSEMBL)

# 5. Maak PWF zoals je al had
pwf <- nullp(
  DEgenes   = all_genes_clean,
  genome    = "hg38",
  id        = "ENSEMBL",
  bias.data = bias.data_clean
)

# 6. Run GO analyse
GO_results <- goseq(
  pwf,
  gene2cat = gene2cat,
  method = "Hypergeometric",
  use_genes_without_cat = TRUE
)

head(GO_results)

#ggplot visualisatie
library(ggplot2)
library(dplyr)

GO_results %>%
  top_n(10, wt = -over_represented_pvalue) %>%
  mutate(hitsPerc = numDEInCat * 100 / numInCat) %>%
  ggplot(aes(x = hitsPerc,
             y = term,
             colour = over_represented_pvalue,
             size = numDEInCat)) +
  geom_point() +
  expand_limits(x = 0) +
  labs(x = "Hits (%)", y = "GO term", colour = "p value", size = "Count")


# GOTERM[[...]] → haalt alle informatie op die bij die GO-ID hoort
GOTERM[[GO_results$category[1]]]





