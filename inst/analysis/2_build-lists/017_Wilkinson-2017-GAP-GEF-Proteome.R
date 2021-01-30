#!/usr/bin/env Rscript

## Scrape synaptic GAP/GEF Proteome from Wilkinson et al., 2017

## Parameters:
script <- "017_Wilkinson-2017-GAP-GEF-Proteome"
short_name <- "wilkinson2017gapgef"
paper <- "https://www.ncbi.nlm.nih.gov/pubmed/28706196"

# Urls to paper and data.
urls <- list(
  S1 = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-017-05588-3/MediaObjects/41598_2017_5588_MOESM2_ESM.xls",
  S2 = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-017-05588-3/MediaObjects/41598_2017_5588_MOESM3_ESM.xls",
  S3 = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-017-05588-3/MediaObjects/41598_2017_5588_MOESM4_ESM.xls",
  S4 = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-017-05588-3/MediaObjects/41598_2017_5588_MOESM5_ESM.xls"
)

renv::load(getrd())

# Imports
library(getPPIs)
library(readxl)

# Additional functions.
devtools::load_all()

# Directories.
root <- getrd()
gmtdir <- file.path(root, "gmt")
datadir <- file.path(root, "data")

# Download the data.
message("Downloading data...")
download.files(urls, quiet = TRUE)
myfiles <- basename(unlist(urls))
names(myfiles) <- names(urls)

# S3 has the data we need.
data <- list()
data[["Syngap1"]] <- read_excel(myfiles[["S3"]], "Syngap1 PSD")
data[["Agap2"]] <- read_excel(myfiles[["S3"]], "Agap2 PSD", skip = 1)
data[["Kalrn"]] <- read_excel(myfiles[["S3"]], "Kalirin PSD")

# Remove temporary files.
unlink(myfiles)

# Map mgi ids to entrez and create gmt.
data <- lapply(data, function(x) {
  entrez <- getIDs(x[["MGI ID"]], from = "mgi", to = "entrez", species = "mouse")
  x$entrez <- entrez
  return(x)
})

# Collect genes.
gene_list <- lapply(data, function(x) unique(x$entrez))
gene_list <- lapply(gene_list, function(x) x[!is.na(x)])

# Summary.
knitr::kable(t(sapply(gene_list, length)))

# Write gmt.
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, paper, myfile)

# Save as rda and generate documentation.
documentDataset(myfile, short_name, Rdir = file.path(root, "R"), datadir)
