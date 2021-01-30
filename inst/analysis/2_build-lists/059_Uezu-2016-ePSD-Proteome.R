#!/usr/bin/env Rscript

# Scrape iPSD proteome from Uezu et al., 2016

# Load renv.
here <- getwd()
root <- dirname(dirname(here))
renv::load(root)

# Urls to data.
data_source <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5432043/"
data_url <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5432043/bin/NIHMS855807-supplement-Table_S3.xlsx"
script <- "059_Uezu-2016-ePSD-Proteome"
short_name <- "ePSD"

# Imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(getPPIs)
  library(data.table)
})

# Functions.
suppressWarnings({
  devtools::load_all()
})

# Directories.
gmtdir <- file.path(root, "gmt")
datadir <- file.path(root, "data")

# Download the data.
myfile <- basename(data_url)
download.file(data_url, myfile, quiet = TRUE)

# Read the data.
data <- read_excel_sheets(myfile, skip = 1)
data <- data[[1]] # There is only one sheet.

# Remove temporary file.
invisible(unlink(myfile))

# Map Uniprot to Entrez.
uniprot <- sapply(strsplit(trimws(data$"UniProt accession"),"-"),"[",1)
entrez <- getIDs(uniprot, from = "uniprot", to = "entrez", species = "mouse")

# Use gene symbols to try mapping unmapped genes.
not_mapped <- is.na(entrez)
symbols <- data$"GeneName"
entrez[not_mapped] <- getIDs(symbols[not_mapped],
  from = "symbol", to = "entrez", species = "mouse"
)

if (sum(is.na(entrez)) != 0) { stop() }

# Add to the data.
data <- tibble::add_column(data, entrez, .after = "UniprotID")

# Status.
nGenes <- length(unique(data$entrez))
message(paste("Collected", nGenes, "ePSD genes."))

# Save as gmt.
gene_list <- list(ePSD = unique(data$entrez))
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, data_source, myfile)

# Save as rda and generate documentation.
documentDataset(myfile, short_name, Rdir = file.path(root, "R"), datadir)
