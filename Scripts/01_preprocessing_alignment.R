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
