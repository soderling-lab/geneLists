#!/usr/bin/env Rscript

## Scrape Phillips et al., 2004 PSD proteome.

## Parameters:
short_name <- "phillips2004psd"
script <- "031_Phillips-2004-PSD-Proteome"
ref_url <- "https://www.ncbi.nlm.nih.gov/pubmed?term=15372503"

# Data from Supplemental table 1 were downloaded as pdf and
# converted to excel. The excel document was processed by
# hand, ans saved as .csv in downloads.

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

# Load the data.
myfile <- file.path(downdir, "Phillips-2004-T1.xlsx")
data <- read_excel(myfile)

# Collect gene identifiers.
ids <- data$Accession
IDs <- file.path(getwd(), "IDs.txt")
writeLines(ids, IDs)

# Call python script to map identifiers.
from <- "REFSEQ_NT_ID"
to <- "P_ENTREZGENEID"
pyScript <- file.path(root, "Py", "mapIds.py")
cmd <- paste(pyScript, IDs, from, to, "--nsteps 2")
result <- system(cmd, intern = TRUE)

# Remove temporary file.
unlink(IDs)

# Parse the result.
entrez <- strsplit(gsub("\\[|\\]|'", "", result), ", ")[[1]]
entrez[entrez == "None"] <- NA

# For missing ids, remove last character.
is_missing <- is.na(entrez)
ids[is_missing] <- str_sub(ids[is_missing], end = -2)

# Try again.
writeLines(ids, IDs)
from <- "REFSEQ_NT_ID"
to <- "P_ENTREZGENEID"
pyScript <- file.path(root, "Py", "mapIds.py")
cmd <- paste(pyScript, IDs, from, to, "--nsteps 2")
result <- system(cmd, intern = TRUE)
unlink(IDs)
entrez <- strsplit(gsub("\\[|\\]|'", "", result), ", ")[[1]]
entrez[entrez == "None"] <- NA

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

# Collect list of genes.
gene_list <- list()
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
