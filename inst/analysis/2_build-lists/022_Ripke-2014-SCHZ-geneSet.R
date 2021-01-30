#!/usr/bin/env Rscript

## Scrape Schizophrenia genes from Fromer et al., 2014
# Table S3 was downloaded and converted to an excel document
# using an online converter:
# https://www.ilovepdf.com/pdf_to_excel

## Parameters:
short_name <- "ripke2014SCHZ"
script <- "022_Ripke-2014-SCHZ-geneSet"
ref_url <- "https://www.ncbi.nlm.nih.gov/pubmed/25056061"

# Load renv:
renv::load(getrd())

# Imports.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(getPPIs)
  library(readxl)
})

# Functions.
devtools::load_all()

# Directories.
root <- getrd()
gmtdir <- file.path(root, "datasets")
datadir <- file.path(root, "data")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "downloads")

# Load the data.
myfile <- file.path(downdir, "ripke-2014-s3.xlsx")
data <- read_excel(myfile)

# Split Protein coding genes column.
# "*" indicate loci where the region does not contain a gene.
# Therefore the nearest gene was listed.
data <- data %>% tidyr::separate_rows("Protein coding genes", sep = " ")
data <- data %>% tidyr::separate_rows("Protein coding genes", sep = "\r\n")

# Drop "*" and NA.
idx <- data$"Protein coding genes" == "*" | is.na(data$"Protein coding genes")
data <- data[!idx, ]

# Map human genes to entrez.
genes <- data$"Protein coding genes"
entrez <- getIDs(genes, from = "symbol", to = "entrez", species = "human")
data <- tibble::add_column(data, entrez, .after = "Protein coding genes")

# Map to mouse homologs.
msEntrez <- getHomologs(entrez, species = "mouse")
data <- tibble::add_column(data, msEntrez, .after = "entrez")

# Drop NA.
data <- data[!is.na(data$msEntrez), ]

# Status
genes <- unique(data$msEntrez)
nGenes <- length(genes)
message(paste("Compiled", nGenes, "mouse genes associated with schizophrenia."))

# Save as gene list.
myfile <- file.path(gmtdir, script)
gene_list <- list(short_name = genes)
write_gmt(gene_list, ref_url, myfile)

# Save as rda and generate documentation.
documentDataset(myfile, short_name, Rdir = file.path(root, "R"), datadir)
