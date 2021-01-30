#!/usr/bin/env Rscript

## Scrape Collins PSD proteome.

## Parameters:
short_name <- "collins2006PSD"
script <- "027_Collins-2006-PSD-Proteome"
ref_url <- "https://www.ncbi.nlm.nih.gov/pubmed?term=16635246"

# Data from Supplemental table 3 were downloaded as pdf and
# converted to excel. The excel document was processed by
# hand and stored in downloads.

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
myfile <- file.path(downdir, "Collins-2006-S3.xlsx")
data <- read_excel(myfile)

# Map uniprot ids with uniprot API.
uniprot <- data$UniProt
IDs <- file.path(getwd(), "IDs.txt")
writeLines(uniprot, "IDs.txt")

# Call python script to map identifiers.
pyScript <- file.path(root, "Py", "mapIds.py")
cmd <- paste(pyScript, IDs, "ACC", "P_ENTREZGENEID", "--nsteps 1")
result <- system(cmd, intern = TRUE)

# Remove temporary file.
invisible(unlink(IDs))

# Parse the result.
entrez <- strsplit(gsub("\\[|\\]|'", "", result), ", ")[[1]]
entrez[entrez == "None"] <- NA

# Try to map missing ids with unigene.
is_missing <- is.na(entrez)
sum(is_missing)
symbol <- data$"Gene name"[is_missing]
entrez[is_missing] <- getIDs(symbol,
  from = "symbol",
  to = "entrez", species = "rat"
)

# Check how many are missing.
is_missing <- is.na(entrez)
n_missing <- sum(is_missing)
message(paste(
  "Percent of genes mapped to stable entrez IDs:",
  round(100 * (1 - (n_missing / length(is_missing))), 3)
))

# Map rat entrez to mouse.
msEntrez <- getHomologs(entrez, species = "mouse")

# Add to data.
data <- tibble::add_column(data, entrez, .after = "Gene name")
data <- tibble::add_column(data, msEntrez, .after = "entrez")

# Drop NA.
data <- data[!is.na(data$msEntrez), ]

# Collect groups of genes.
data_list <- split(data, data$"Protein Subclass")

# Collect list of genes.
gene_list <- lapply(data_list, function(x) unique(x$msEntrez))
gene_list[["All"]] <- unique(data$msEntrez)

# Summary
df <- data.frame("Proteins" = sapply(gene_list, length))
df <- tibble::add_column(df, "Function" = rownames(df), .before = 1)
knitr::kable(df, row.names = FALSE)

# Save as gene list.
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, ref_url, myfile)

# Save as rda and generate documentation.
documentDataset(myfile, short_name, Rdir = file.path(root, "R"), datadir)
