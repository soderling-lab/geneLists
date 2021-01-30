#!/usr/bin/env Rscript

## Building a ASD-risk gene dataset from the WES data from
# Sanders et al., 2015.
# Supplemental data were downloaded from:
# https://www.sciencedirect.com/science/article/pii/S0896627315007734?via%3Dihub#app2

renv::load(getrd())

## User parameters:
script <- "007_Sanders-2015-ASD-geneSet"
short_name <- "sanders2015ASD"
gene_source <- "Sanders_et_al_2015"
data_source <- "https://www.ncbi.nlm.nih.gov/pubmed/26402605"
data_description <- "ASD-associated genes."

# Imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(getPPIs)
  library(data.table)
})

# Load functions.
devtools::load_all()

# Directories.
root <- getrd()
gmtdir <- file.path(root, "gmt")
funcdir <- file.path(root, "R")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "downloads")

# Download the data.
url <- "https://ars.els-cdn.com/content/image/1-s2.0-S0896627315007734-mmc7.xlsx"
myfile <- file.path(downdir, basename(url))
download.file(url, destfile = myfile, quiet = TRUE)

# Load the data.
data <- read_excel(myfile)

# Keep genes with FDR < 0.1 from these columns.
cols <- c(
  "65genes_tadaFdrAscSscExomeSscAgpSmallDel",
  "59genes_tadaFdrAscSscExome"
)
keep <- data[cols[1]] != "0" | data[cols[2]] != "0"
data <- subset(data, keep)

# Map gene symbols to entrez.
symbols <- data$RefSeqGeneName
entrez <- getIDs(symbols, from = "symbol", to = "entrez", species = "human")

# Map missing by hand.
not_mapped <- which(is.na(entrez))
entrez[not_mapped] <- c(51111, 55914)

# Add to data.
data <- tibble::add_column(data, "entrez" = entrez, .after = "TadaGeneName")

# Map human genes in to their mouse homologs.
msEntrez <- getHomologs(entrez, taxid = 10090)
data <- tibble::add_column(data, msEntrez = msEntrez, .after = "entrez")

# Status.
n_mapped <- sum(!is.na(msEntrez))
n_genes <- length(entrez)
message(paste(
  n_mapped, "of", n_genes, "human ASD-risk",
  "genes were successfully mapped to mouse homologs."
))

# Remove rows with unmapped genes.
data <- subset(data, !is.na(msEntrez))

# Save table.
myfile <- file.path(tabsdir, paste0(script, ".csv"))
fwrite(data, myfile)

# Save as gmt.
gmt_list <- list("ASD" = data$msEntrez)
gmt_file <- file.path(gmtdir, script)
gmt_source <- data_source
write_gmt(gmt_list, gmt_source, gmt_file)

# Save as rda and generate documentation.
documentDataset(gmt_file, short_name, Rdir = file.path(root, "R"), datadir)
