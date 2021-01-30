#!/usr/bin/env Rscript

## ARGS:
script <- "054_Spence-2019-Wrp-Proteome"
short_name <- "spence2019"
tags <- c("BioID", "proteomics", "Wrp", "nascent synapse")
ref_url <- "https://pubmed.ncbi.nlm.nih.gov/30674877/"
data_url <- "downloads/Wrp-BirA-enriched.xlsx"

# Load renv.
here <- getwd()
root <- dirname(dirname(here))
renv::load(root)

# Imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(getPPIs)
  library(TBmiscr)
  library(data.table)
})

# Load functions in root/R.
TBmiscr::load_all()

# Directories.
root <- getrd()
gmtdir <- file.path(root, "gmt")
datadir <- file.path(root, "data")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "downloads")

# Load the data.
message("\nLoading the data from file.")
myfile <- file.path(root, data_url)
wrp_data <- read_excel(myfile)

# Try mapping genes to entrez.
symbols <- wrp_data$Label
entrez <- getPPIs::getIDs(symbols, from = "symbol", to = "entrez", species = "mouse")

# Map missing ids by hand.
message("\nMapping missing ids by hand using MGI website.")
missing <- is.na(entrez)
mapped_by_hand <- c(
  "Azi1" = 12009,
  "Bai3" = 210933,
  "Erbb2ip" = 59079,
  "Kiaa1549" = 330286,
  "Lphn1" = 330814,
  "Lphn2" = 99633,
  "Lphn3" = 319387,
  "Lppr2" = 235044,
  "Lppr4" = 229791,
  "Lrrc16b" = 268747,
  "Odz1" = 23963,
  "Odz2" = 23964,
  "Odz3" = 23965,
  "Odz4" = 23966,
  "Pvrl3" = 58998
)
entrez[missing] <- mapped_by_hand[names(entrez[missing])]

if (any(is.na(entrez))) {
  stop("Trouble mapping gene names to entrez ids.")
}
if (length(unique(entrez)) != length(entrez)) {
  stop("duplicates!")
}

# Status.
ngenes <- length(entrez)
message(paste("\nCompiled", ngenes, "mouse proteins in the WRP P5 iBioID proteome."))

# Save this as a gene list.
gene_list <- list("wrp_P5-interactome" = entrez)

# Save as gmt, and then save as rda and generate documentation.
Rdir <- file.path(root, "R")
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, ref_url, myfile)
documentDataset(myfile, short_name, Rdir, datadir)
