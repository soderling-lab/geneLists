#!/usr/bin/env Rscript

## Scrape Lelievald et al., 2016 ID genes.

## Parameters:
short_name <- "lelieveld2016ID"
script <- "038_Lelieveld-2016-ID-geneSet"
ref_url <- "https://www.ncbi.nlm.nih.gov/pubmed/27479843"
data_url <- c(
  "https://static-content.springer.com/esm/art%3A10.1038%2Fnn.4352/MediaObjects/41593_2016_BFnn4352_MOESM23_ESM.xlsx",
  "https://static-content.springer.com/esm/art%3A10.1038%2Fnn.4352/MediaObjects/41593_2016_BFnn4352_MOESM25_ESM.xls"
)
pmid <- "27479843"

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
download.files(data_url, destdir = downdir, quiet = TRUE)
myfiles <- file.path(downdir, destdir = basename(data_url))
data <- list()
data[["ID Genes"]] <- read_excel(myfiles[1])
data[["DNM Enriched"]] <- read_excel(myfiles[2])

# Collect Diagnostic genes.
df <- data[[1]]
entrez <- getIDs(df$Genes,
  from = "symbol", to = "entrez", species = "human"
)
msEntrez <- getHomologs(entrez, species = "mouse")
idx <- df$Diagnostics

# Compile list of genes.
gene_list <- list()
gene_list[["ID Diagnostic Genes"]] <- unique(msEntrez[idx])

# Get top DNM genes.
df <- data[[2]]
entrez <- getIDs(df$gene.name,
  from = "symbol", to = "entrez", species = "human"
)
msEntrez <- getHomologs(entrez, species = "mouse")
idx <- df$significant.based.on.FDR < 0.1
gene_list[["DNM Enriched"]] <- unique(msEntrez[idx])
gene_list[["Combined"]] <- unique(unlist(gene_list))

# Summary
df <- data.frame("N Genes" = sapply(gene_list, length))
df <- tibble::add_column(df, "Class" = rownames(df), .before = 1)
knitr::kable(df, row.names = FALSE)

# Save as gene list.
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, ref_url, myfile)

# Save as rda and generate documentation.
documentDataset(myfile, short_name, Rdir = file.path(root, "R"), datadir)
