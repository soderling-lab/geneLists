#!/usr/bin/env Rscript

## Scrape Lee et al., 2017 Shank3 Zinc-induced interactome.

## Parameters:
short_name <- "lee2017shank3zinc"
script <- "039_Lee-2017-Shank3-Zn-Proteome"
ref_url <- "https://www.ncbi.nlm.nih.gov/pubmed/29111324/"
data_url <- "https://ars.els-cdn.com/content/image/1-s2.0-S0006291X17321289-mmc1.xlsx"
pmid <- "29111324"

# Get data from sheet 2, the authors zinc induced interactome.

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
download.file(data_url, myfile, quiet = TRUE)
data <- read_excel(myfile, sheet = 2, skip = 1)

# Split Accession rows with multiple entries.
data <- tidyr::separate_rows(data, Accession, sep = ";")
data <- data %>% filter(Accession != "")

# Map (mouse) uniprot accession ids to entrez.
uniprot <- data$Accession
entrez <- getIDs(uniprot, from = "uniprot", to = "entrez", species = "mouse")
is_missing <- is.na(entrez)
entrez[is_missing] <- getIDs(data$Mouse[is_missing],
  from = "symbol", to = "entrez", species = "mouse"
)

# Add to data.
data <- tibble::add_column(data, entrez, .after = "Accession")

# Compile list of genes.
gene_list <- list()
gene_list[["Zn-Induced Interactome"]] <- unique(data$entrez)

# Summary
df <- data.frame("N Genes" = sapply(gene_list, length))
df <- tibble::add_column(df, "Class" = rownames(df), .before = 1)
knitr::kable(df, row.names = FALSE)

# Save as gene list.
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, ref_url, myfile)

# Save as rda and generate documentation.
documentDataset(myfile, short_name, Rdir = file.path(root, "R"), datadir)
