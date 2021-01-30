#!/usr/bin/env Rscript

## Scrape Collins et al., 2005 Synaptic phospho-proteome.

## Parameters:
short_name <- "collins2005phospho"
script <- "028_Collins-2005-Synaptic-Phospho-Proteome"
ref_url <- "https://www.ncbi.nlm.nih.gov/pubmed?term=15572359"

# Data from Supplemental table 1 were downloaded as pdf and
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
myfile <- file.path(downdir, "Collins-2005-S1.xlsx")
data <- read_excel(myfile)

# Cleanup the columns.
colnames(data) <- gsub("\r|\n", "", colnames(data))

# Load kinase annotations.
myfile <- file.path(downdir, "Collins-2005-Kinases.csv")
kinase_map <- fread(myfile)
kinases <- kinase_map$Kinase
names(kinases) <- kinase_map$Abbreviation

# Separate rows with multiple putative kinases.
data <- tidyr::separate_rows(data, "Known cognate kinases (literature)", sep = ", ")
data <- tidyr::separate_rows(data, "Scansite phosphorylation hits", sep = ", ")

# Get predicted kinase.
data$"Predicted Kinase" <- kinases[data$"Scansite phosphorylation hits"]
data$"Scansite phosphorylation hits" <- NULL

# Map uniprot ids with uniprot API.
uniprot <- data$Accession
IDs <- file.path(getwd(), "IDs.txt")
writeLines(uniprot, "IDs.txt")

# Call python script to map identifiers.
pyScript <- file.path(root, "Py", "mapIds.py")
cmd <- paste(pyScript, IDs, "ACC", "P_ENTREZGENEID", "--nsteps 1")
result <- system(cmd, intern = TRUE)

# Remove temporary file.
unlink(IDs)

# Parse the result.
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

# Convert NA to unknown.
data$"Predicted Kinase"[is.na(data$"Predicted Kinase")] <- "Unknown"
data$"Known cognate kinases (literature)"[is.na(data$"Known cognate kinases (literature)")] <- "Unknown"

# Group by predicted kinase.
list1 <- split(data, data$"Predicted Kinase")
x <- list1$Unknown
list1[["Unknown"]] <- NULL
names(list1) <- paste("Predicted", names(list1), "substrates")

# Group by known kinase.
list2 <- split(data, data$"Known cognate kinases (literature)")
y <- list2$Unknown
list2[["Unknown"]] <- NULL
names(list2) <- paste("Known", names(list2), "substrates")

# Combine and add unknown substrates.
data_list <- c(list1, list2)
data_list[["Unknown substrates"]] <- rbind(list1$Unknown, list2$Unknown)

# Collect list of genes.
gene_list <- lapply(data_list, function(x) unique(x$msEntrez))
gene_list[["All Phosphoproteins"]] <- unique(data$msEntrez)

# Summary
df <- data.frame("Proteins" = sapply(gene_list, length))
df <- tibble::add_column(df, "Substrate Class" = rownames(df), .before = 1)
knitr::kable(df, row.names = FALSE)

# Save as gene list.
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, ref_url, myfile)

# Save as rda and generate documentation.
documentDataset(myfile, short_name, Rdir = file.path(root, "R"), datadir)
