#!/usr/bin/env Rscript

## Parameters:
short_name <- "trinidad2005phosphoPSD"
script <- "033_Trinidad-2005-Phospho-PSD-Proteome"
ref_url <- "http://www.ncbi.nlm.nih.gov/pubmed?term=15748150"

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
myfile <- file.path(downdir, "Trinidad-2005-S1.xlsx")
data <- read_excel(myfile)

# Collect ids.
ids <- data$Accession
entrez <- getIDs(ids, from = "refseq", to = "entrez", species = "mouse")

# Check how many are missing.
is_missing <- is.na(entrez)
n_missing <- sum(is_missing)
message(paste(
  "Percent of genes mapped to stable entrez IDs:",
  round(100 * (1 - (n_missing / length(is_missing))), 3)
))

# Add to data.
data <- tibble::add_column(data, entrez, .after = "Accession")

# Drop NA.
data <- data[!is.na(data$entrez), ]

# Group genes into functional classes.
data_list <- split(data, data$"Protein Function")

# Collect list of genes.
gene_list <- lapply(data_list, function(x) unique(x$entrez))

# Summary
df <- data.frame("Proteins" = sapply(gene_list, length))
df <- tibble::add_column(df, "Protein Class" = rownames(df), .before = 1)
knitr::kable(df, row.names = FALSE)

# Save as gene list.
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, ref_url, myfile)

# Save as rda and generate documentation.
documentDataset(myfile, short_name, Rdir = file.path(root, "R"), datadir)
