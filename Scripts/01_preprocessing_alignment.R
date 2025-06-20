# Uitpakken en indexeren (Scripts/1_preprocessing_indexing.R)
setwd("C:/Users/Elda Mehmedagic/School")

# Uitpakken van ZIP-bestand
install.packages("R.utils")  # Alleen bij eerste keer nodig
library(R.utils)
unzip("Data_RA_raw.zip", exdir = "Data_RA_raw")

# Indexeren van het referentiegenoom
install.packages("BiocManager")
BiocManager::install("Rsubread")
library(Rsubread)

buildindex( basename = "ref_human",
            reference = "GCF_000001405.40_GRCh38.p14_genomic.fna",
            memory = 7000,
            indexSplit = TRUE)

# FASTQ-bestanden alignen naar BAM
align.ra1 <- align(index = "ref_human", readfile1 = "SRR4785819_1_subset40k.fastq", readfile2 = "SRR4785819_2_subset40k.fastq", output_file = "ra1.BAM")
align.ra2 <- align(index = "ref_human", readfile1 = "SRR4785820_1_subset40k.fastq", readfile2 = "SRR4785820_2_subset40k.fastq", output_file = "ra2.BAM")
align.ra3 <- align(index = "ref_human", readfile1 = "SRR4785828_1_subset40k.fastq", readfile2 = "SRR4785828_2_subset40k.fastq", output_file = "ra3.BAM")
align.ra4 <- align(index = "ref_human", readfile1 = "SRR4785831_1_subset40k.fastq", readfile2 = "SRR4785831_2_subset40k.fastq", output_file = "ra4.BAM")
align.ra5 <- align(index = "ref_human", readfile1 = "SRR4785979_1_subset40k.fastq", readfile2 = "SRR4785979_2_subset40k.fastq", output_file = "ra5.BAM")
align.ra6 <- align(index = "ref_human", readfile1 = "SRR4785980_1_subset40k.fastq", readfile2 = "SRR4785980_2_subset40k.fastq", output_file = "ra6.BAM")
align.ra7 <- align(index = "ref_human", readfile1 = "SRR4785986_1_subset40k.fastq", readfile2 = "SRR4785986_2_subset40k.fastq", output_file = "ra7.BAM")
align.ra8 <- align(index = "ref_human", readfile1 = "SRR4785988_1_subset40k.fastq", readfile2 = "SRR4785988_2_subset40k.fastq", output_file = "ra8.BAM")