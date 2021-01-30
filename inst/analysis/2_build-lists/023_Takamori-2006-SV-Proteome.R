#!/usr/bin/env Rscript

## Scrape Synaptic vesicle proteome from Takamori et al., 2006

## Parameters:
short_name <- "takamori2006SV"
script <- "023_Takamori-2006-SV-Proteome"
ref_url <- "https://www.ncbi.nlm.nih.gov/pubmed/17110340"

# Load renv.
here <- getwd()
root <- dirname(dirname(here))
renv::load(root)

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
myfile <- file.path(downdir, "Takamori-2006-S1.xlsx")
data <- read_excel(myfile)

# Map refseq.
entrez <- getIDs(data$RefSeq, from = "refseq", to = "entrez", species = "rat")
is_missing <- is.na(entrez)

# Use gi to map missing ids to entrez.
gi <- sapply(strsplit(data$"Accession No.", "gi\\|"), "[", 2)
IDs <- file.path(getwd(), "IDs.txt")
writeLines(gi[is_missing], IDs)

# Call python script to map identifiers.
pyScript <- file.path(root, "Py", "mapIds.py")
cmd <- paste(pyScript, IDs, "P_GI", "P_ENTREZGENEID")
result <- system(cmd, intern = TRUE)

# Remove temporary file.
unlink(IDs)

# Parse the result.
entrez[is_missing] <- strsplit(gsub("\\[|\\]|'", "", result), ", ")[[1]]

# Replace None with NA and coerce to numeric.
entrez[entrez == "None"] <- NA

# Status.
not_mapped <- 100 * sum(is.na(entrez)) / length(entrez)
message(paste("Percent unmapped genes:", round(not_mapped, 3)))

# Map rat entrez to mouse.
msEntrez <- getHomologs(entrez, species = "mouse")

# Add to data.
data$entrez <- entrez
data$msEntrez <- msEntrez

# Drop NA.
data <- data[!is.na(data$msEntrez), ]

# Status
genes <- unique(data$msEntrez)
nGenes <- length(genes)
message(paste("Compiled", nGenes, "mouse SV-associated proteins."))

# Group into gene groups.
data_list <- split(data, data$Class)

# Collect lists of genes.
gene_list <- lapply(data_list, function(x) unique(x$msEntrez))
gene_list[["All"]] <- genes

# Summary
df <- data.frame("N Genes" = sapply(gene_list, length))
df <- tibble::add_column(df, "Class" = rownames(df), .before = 1)
knitr::kable(df, row.names = FALSE)

# Save as gene list.
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, ref_url, myfile)

# Save as rda and generate documentation.
documentDataset(myfile, short_name, Rdir = file.path(root, "R"), datadir)
