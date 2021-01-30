#!/usr/bin/env Rscript

## Scrape Jordan et al., 2004 PSD proteome.

## Parameters:
short_name <- "jordan2004psd"
script <- "029_Jordan-2004-PSD-Proteome"
ref_url <- "http://www.ncbi.nlm.nih.gov/pubmed?term=15169875"

renv::load(getrd())

# Data from Supplemental table 1 were downloaded as pdf and
# converted to excel. The excel document was processed by
# hand, ans saved as .csv in downloads.

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
myfile <- file.path(downdir, "Jordan-2004-T1.csv")
data <- fread(myfile)

# Map refseq ids.
refseq <- data$Accession

# Map refseq -- mixture of rat and mouse.
entrez <- getIDs(refseq, from = "refseq", to = "entrez", species = "rat")
entrez[is.na(entrez)] <- getIDs(refseq[is.na(entrez)],
  from = "refseq", to = "entrez", species = "mouse"
)
names(entrez) <- refseq

# Map uniprot.
idx <- grepl("NP_|XP_", data$Accession)
entrez[!idx] <- getIDs(refseq[!idx], from = "uniprot", to = "entrez", species = "mouse")


# Check how many are missing.
is_missing <- is.na(entrez)
n_missing <- sum(is_missing)
message(paste(
  "Percent of genes mapped to stable entrez IDs:",
  round(100 * (1 - (n_missing / length(is_missing))), 3)
))

# Map entrez ids to mouse homologs.
msEntrez <- getHomologs(entrez, species = "mouse")

# Add to data.
data <- tibble::add_column(data, entrez, .after = "Accession")
data <- tibble::add_column(data, msEntrez, .after = "entrez")

# Drop NA.
data <- data[!is.na(data$msEntrez), ]

# Group by protein class/function.
data_list <- split(data, data$"Protein Class")

# Collect list of genes.
gene_list <- lapply(data_list, function(x) unique(x$msEntrez))
gene_list[["All"]] <- unique(data$msEntrez)

# Summary
df <- data.frame("Proteins" = sapply(gene_list, length))
df <- tibble::add_column(df, "Protein Class" = rownames(df), .before = 1)
knitr::kable(df, row.names = FALSE)

# Save as gene list.
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, ref_url, myfile)

# Save as rda and generate documentation.
documentDataset(myfile, short_name, Rdir = file.path(root, "R"), datadir)
