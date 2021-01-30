#!/usr/bin/env Rscript

## Building a ASD-risk gene dataset from the WES data from Satterstrom et al.,
# 2020. Supplemental data were downloaded from:
# https://www.cell.com/cell/fulltext/S0092-8674(19)31398-4?rss=yes#secsectitle0385

renv::load(getrd())

# Input args:
script <- "008_Satterstrom-2020-ASD-geneSet"
short_name <- "satterstrom2020ASD"
gene_source <- "Satterstrom_et_al_2020"
data_source <- "https://www.biorxiv.org/content/10.1101/484113v3"
output_file <- "mouse_Satterstrom_ASD_geneSet.RData"
description <- "ASD/DDID-associated genes."

# Functions.
devtools::load_all()

# Imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(getPPIs)
  library(data.table)
})

# Directories.
root <- getrd()
gmtdir <- file.path(root, "gmt")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "downloads")

# Download the Satterstrom data.
# We need the data from S2.
url <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867419313984-mmc2.xlsx"
myfile <- file.path(downdir, basename(url))
download.file(url, destfile = myfile, quiet = TRUE)

# Read the data.
sheet_names <- excel_sheets(myfile)
all_data <- read_excel_sheets(myfile)

# Get ASD genes.
data <- all_data[["102_ASD"]]

# Remove the last three rows--these just contain some comments.
data <- data[!is.na(data$entrez_id), ]

# Get human entrez ids.
hsEntrez <- data$entrez_id

# Status.
n_genes <- sum(!is.na(hsEntrez))
check <- n_genes == 102
if (!check) {
  stop("Warning, problem parsing excel data.")
}

# Map human genes in to their mouse homologs.
msEntrez <- getHomologs(hsEntrez, taxid = 10090)
data <- tibble::add_column(data, msEntrez = msEntrez, .after = "entrez_id")

# Status.
n_mapped <- sum(!is.na(data$msEntrez))
message(paste(
  n_mapped, "of", n_genes, "human ASD-risk",
  "genes were successfully mapped to mouse homologs."
))

# Remove rows with unmapped genes.
data <- subset(data, !is.na(msEntrez))

# Save table.
myfile <- file.path(tabsdir, paste0(script, ".csv"))
fwrite(data, myfile)

# Write as gmt.
gmt_list <- list("ASD" = data$msEntrez)
gmt_source <- "https://www.biorxiv.org/content/10.1101/484113v3"
gmt_file <- file.path(gmtdir, script)
write_gmt(gmt_list, gmt_source, gmt_file)

# Save as rda and generate documentation.
documentDataset(gmt_file, short_name, Rdir = file.path(root, "R"), datadir)
