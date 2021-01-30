#!/usr/bin/env Rscript

## Scrape ASD genes from Iossifov et al., 2014

## Parameters:
script <- "018_Iossifov-2014-ASD-geneSet"
short_name <- "iossifov2014ASD"
reference_url <- "https://www.ncbi.nlm.nih.gov/pubmed/25363768"
data_url <- "https://static-content.springer.com/esm/art%3A10.1038%2Fnature13908/MediaObjects/41586_2014_BFnature13908_MOESM117_ESM.zip"

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
zipfile <- file.path(downdir, basename(data_url))
download.file(data_url, destfile = zipfile, quiet = TRUE)

# Unzip.
unzip(zipfile, exdir = downdir)

# List all zipped files.
myfiles <- file.path(downdir, unzip(zipfile, list = TRUE)$Name)
myfiles <- myfiles[grepl("xlsx", myfiles)]

# Remove zipped file.
unlink(zipfile)

# Load all the data.
alldat <- lapply(myfiles, read_excel_sheets)
alldat <- unlist(alldat, recursive = FALSE)

# The authors focus on De novo likely gene-disrupting mutations.
data <- alldat[["ggg_clean"]]
data <- data %>% filter(dnv_LGDs_prb != 0) # De novo LGD mutations in affected individuals (probands)

# Map human genes to entrez ids.
entrez <- getIDs(data$gene, from = "symbol", to = "entrez", species = "human")
data <- tibble::add_column(data, entrez, .after = "gene")

# Remove NA.
data <- data[!is.na(data$entrez), ]

# Map human genes to mouse homologs.
msEntrez <- getHomologs(data$entrez, species = "mouse")
data <- tibble::add_column(data, msEntrez, .after = "entrez")

# Remove NA.
data <- data[!is.na(data$msEntrez), ]

# Save table.
myfile <- file.path(tabsdir, paste0(script, ".csv"))
fwrite(data, myfile)

# status
genes <- unique(data$msEntrez)
nGenes <- length(genes)
message(paste(
  "Compiled", nGenes,
  "genes with likely gene-disrupting (LGD) mutations in idividuals with ASD."
))

# Save as gene list.
myfile <- file.path(gmtdir, script)
gene_list <- list("ASD" = genes)
write_gmt(gene_list, reference_url, myfile)

# Save as rda and generate documentation.
documentDataset(myfile, short_name, Rdir = file.path(root, "R"), datadir)
