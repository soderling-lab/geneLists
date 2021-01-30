#!/usr/bin/env Rscript

## Parameters:
short_name <- "chen2014AMPARs"
script <- "040_Chen-2014-AMPAR-Proteome"
ref_url <- "https://www.ncbi.nlm.nih.gov/pubmed/25337787"
data_url <- "https://pubs.acs.org/doi/suppl/10.1021/pr500697b/suppl_file/pr500697b_si_003.xlsx"
pmid <- "2533787"

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
myfile <- file.path(downdir, basename(data_url))
download.file(data_url, destfile = myfile)
data <- read_excel(myfile, skip = 6)

# Remove rows in which protein was not identified by any of the IPs.
out <- data %>%
  dplyr::select(c("Cerebellum", "Cortex", "Hippocampus")) %>%
  apply(1, function(x) all(is.na(x)))
data <- data[!out, ]

# Separate rows with multiple uniprot accessions or gene names.
data <- tidyr::separate_rows(data, "Entry name (UniProtKB/TrEMBL)", sep = "; ")
data <- tidyr::separate_rows(data, "Gene name", sep = ";")

# Collect uniprot ids.
uniprot <- sapply(strsplit(data$"Entry name (UniProtKB/TrEMBL)", "\\|"), "[", 2)
uniprot <- gsub("*-[0-9]{1}$", "", uniprot)
data <- tibble::add_column(data, uniprot, .after = "Gene name")

# Map mouse uniprot to entrez.
entrez <- getIDs(uniprot, from = "uniprot", to = "entrez", species = "mouse")
is_missing <- is.na(entrez)
sum(is_missing)

# For missing ids, map human symbols to mouse entrez.
symbols <- data$"Gene name"[is_missing]
hsEntrez <- getIDs(symbols, from = "symbol", to = "entrez", species = "human")
entrez[is_missing] <- getHomologs(hsEntrez, species = "mouse")
is_missing <- is.na(entrez)
sum(is_missing)

# Add to data.
data <- tibble::add_column(data, entrez, .after = "uniprot")

# Drop na.
data <- data[!is.na(data$entrez), ]

# Compile list of genes.
gene_list <- list()
gene_list[["AMPAR Proteome"]] <- unique(data$entrez)

# Summary
df <- data.frame("N Genes" = sapply(gene_list, length))
df <- tibble::add_column(df, "Class" = rownames(df), .before = 1)
knitr::kable(df, row.names = FALSE)

# Save as gene list.
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, ref_url, myfile)

# Save as rda and generate documentation.
documentDataset(myfile, short_name, Rdir = file.path(root, "R"), datadir)
