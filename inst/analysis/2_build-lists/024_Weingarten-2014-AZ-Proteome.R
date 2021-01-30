#!/usr/bin/env Rscript

## Scrape Proteome of presynaptic active zone.
# Due to the authors formating of gene names and lack of
# stable identifiers, we are only able to collect ~50% of the
# proteins they report.

## Parameters:
short_name <- "weingarten2014AZ"
script <- "024_Weingarten-2014-AZ-Proteome"
ref_url <- "https://www.ncbi.nlm.nih.gov/pubmed/24534009"
data_url <- "https://ars.els-cdn.com/content/image/1-s2.0-S1044743114000220-mmc1.xlsx"

renv::load(getrd())

# Imports.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(getPPIs)
  library(readxl)
  library(stringr)
})

# Functions.
devtools::load_all()

# Directories.
root <- getrd()
gmtdir <- file.path(root, "datasets")
datadir <- file.path(root, "data")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "downloads")

# Download the data.
myfile <- file.path(downdir, basename(data_url))
download.file(data_url, myfile, quiet = TRUE)
data <- read_excel(myfile)

# Data is from mouse!
# Convert to title case.
genes <- str_to_title(data$"UniProt Acc.")
# Limited success with other approaches.
# genes <- data$"UniProt Acc."
# entrez <- mapUniprot(genes,from="GENENAME",to="P_ENTREZGENEID",root)
# entrez <- mapUniprot(genes,from="GENENAME",to="P_ENTREZGENEID",root)
entrez <- getIDs(genes, from = "symbol", to = "entrez", species = "mouse")

# Status.
not_mapped <- 100 * sum(is.na(entrez)) / length(entrez)
message(paste("Percent unmapped genes:", round(not_mapped, 3)))

# Add to data.
data$entrez <- entrez

# Drop NA.
data <- data[!is.na(data$entrez), ]

# Status
genes <- unique(data$entrez)
nGenes <- length(genes)
message(paste("Compiled", nGenes, "mouse presynaptic AZ proteins."))

# Save as gene list.
gene_list <- list()
gene_list[[short_name]] <- genes
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, ref_url, myfile)

# Save as rda and generate documentation.
documentDataset(myfile, short_name, Rdir = file.path(root, "R"), datadir)
