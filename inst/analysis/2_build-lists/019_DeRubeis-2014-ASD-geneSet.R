#!/usr/bin/env Rscript

## Scrape ASD genes from Iossifov et al., 2014

## Parameters:
short_name <- "derubeis2014ASD"
script <- "019_DeRubeis-2014-ASD-geneSet"
ref_url <- "https://www.ncbi.nlm.nih.gov/pubmed/25363760"
data_url <- "https://static-content.springer.com/esm/art%3A10.1038%2Fnature13772/MediaObjects/41586_2014_BFnature13772_MOESM40_ESM.xlsx"

# Load renv:
renv::load(getrd())

# Imports.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(getPPIs)
})

# Functions.
devtools::load_all()

# Directories.
root <- getrd()
gmtdir <- file.path(root, "gmt")
datadir <- file.path(root, "data")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "downloads")

# Download the data.
myfile <- file.path(downdir, basename(data_url))
download.file(data_url, destfile = myfile, quiet = TRUE)

# Load the data.
data <- read_excel_sheets(myfile)
data <- data[[1]]

# The authors focus on genes with an FDR < 0.3.
data <- data %>% filter(qvalue < 0.3)

# Map human genes to entrez.
entrez <- getIDs(data$Gene, from = "symbol", to = "entrez", species = "human")
data <- tibble::add_column(data, entrez, .after = "Gene")

# Map to mouse.
msEntrez <- getHomologs(data$entrez, species = "mouse")
data <- tibble::add_column(data, msEntrez, .after = "entrez")

# Remove NA.
data <- data[!is.na(data$msEntrez), ]

# Save table.
myfile <- file.path(tabsdir, paste0(script, ".csv"))
fwrite(data, myfile)

# Status
genes <- unique(data$msEntrez)
nGenes <- length(genes)
message(paste("Compiled", nGenes, "genes associated with ASD."))

# Save as gene list.
myfile <- file.path(gmtdir, script)
gene_list <- list("ASD" = genes)
write_gmt(gene_list, ref_url, myfile)

# Save as rda and generate documentation.
documentDataset(myfile, short_name, Rdir = file.path(root, "R"), datadir)
