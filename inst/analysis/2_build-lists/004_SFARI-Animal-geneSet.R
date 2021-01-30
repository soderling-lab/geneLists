#!/usr/bin/env Rscript

## Building a SFARI disease-associated gene collection.

renv::load(getrd())

# Which dataset?
script <- "004_SFARI-Animal-geneSet"
dataset <- "Animal" # SFARI or Animal
short_name <- "sfariAnimal"

# Imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(getPPIs)
  library(data.table)
})

# Load functions.
devtools::load_all()

# Directories.
root <- getrd()
gmtdir <- file.path(root, "gmt")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "downloads")

# Load the data.
datasets <- c(
  SFARI = "SFARI-Gene_genes",
  Animal = "SFARI-Gene_animal-genes"
)
myfile <- list.files(downdir,
  pattern = datasets[dataset][1],
  full.names = TRUE
)
data <- data.table::fread(myfile)
colnames(data) <- gsub("-", "_", colnames(data))

# Map human gene symbols to Entrez.
genes <- data$gene_symbol
entrez <- getIDs(genes, from = "symbol", to = "entrez", species = "human")
data <- tibble::add_column(data, "entrez_id" = entrez, .after = 4)
data <- data[entrez_id != "NA"]

# Map human genes in to their mouse homologs.
hsEntrez <- data$entrez_id
nHsGenes <- length(unique(hsEntrez))
msEntrez <- getHomologs(hsEntrez, taxid = 10090)
data <- tibble::add_column(data, msEntrez = msEntrez, .after = 5)
# Remove rows with unmapped genes.
data <- data[msEntrez != "NA"]

# Status report.
nGenes <- length(unique(data$msEntrez))
message(paste0(
  "Compiled ", nGenes, " mouse genes associated with ",
  "Autism spectrum disorders!"
))

# Write data to file.
myfile <- file.path(tabsdir, paste0(script, ".csv"))
data.table::fwrite(data, myfile)

# Write to gmt.
gmt_list <- list("ASD" = data$msEntrez)
gmt_source <- "https://gene.sfari.org/tools/"
gmt_file <- file.path(gmtdir, paste0(script, ".gmt"))
write_gmt(gmt_list, gmt_source, gmt_file)

# Save as rda and generate documentation.
documentDataset(gmt_file, short_name, Rdir = file.path(root, "R"), datadir)
