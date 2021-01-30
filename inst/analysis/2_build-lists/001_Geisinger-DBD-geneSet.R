#!/usr/bin/env Rscript

# Building a developmental brain disorder-associated gene set
# collection.

# Load renv:
renv::load(getrd())

## User parameters:
dataset <- "All"
script <- "001_Geisinger-DBD-geneSet"
short_name <- "geisingerDBD"
tags <- c("DBD", "ASD")

# Imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(getPPIs)
  library(data.table)
})

# Load functions in R/
suppressWarnings({
  devtools::load_all()
})

# Directories.
root <- getrd()
gmtdir <- file.path(root, "gmt")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "downloads")

# Downoad the raw data.
base_url <- "http://dbd.geisingeradmi.org/downloads"
datasets <- c(
  LOF = "Full-Gene-Table-Data",
  Missense = "Full-Missense-Table-Data",
  All = "DBD-Genes-Full-Data"
)
myurl <- file.path(base_url, paste0(datasets[dataset], ".csv"))
myfile <- paste0(datasets[dataset], ".csv")
download.file(myurl, destfile = basename(myurl), quiet = TRUE)
raw_data <- fread(myfile)
unlink(myfile)

# Map human gene symbols to Entrez.
genes <- unique(raw_data$Gene)
entrez <- getIDs(genes, from = "symbol", to = "entrez", species = "human")
names(entrez) <- genes

# Map missing Entrez IDs by hand.
not_mapped <- genes[is.na(entrez)]
mapped_by_manual_search <- c("HIST1H1E" = 3006, "HIST1H4B" = 8366, "KIF1BP" = 26128)
entrez[not_mapped] <- mapped_by_manual_search[names(entrez[not_mapped])]

# Check.
check <- sum(is.na(entrez)) == 0
if (!check) {
  stop()
}

# Add human entrez IDs to data.
data <- tibble::add_column(raw_data,
  "hsEntrez" = entrez[raw_data$Gene],
  .after = "Gene"
)

# Map human genes in to their mouse homologs.
hsEntrez <- data$hsEntrez
msEntrez <- getHomologs(hsEntrez, species = "mouse")
data <- tibble::add_column(data, msEntrez = msEntrez, .after = "hsEntrez")

# Remove rows with unmapped genes.
data <- data[msEntrez != "NA"]

# Add concise disorder association column.
disorders <- c(
  "ID/DD", "Autism", "Epilepsy", "ADHD",
  "Schizophrenia", "Bipolar Disorder"
)
idy <- sapply(disorders, function(x) grep(x, colnames(data)))
anno <- apply(data[, idy, with = FALSE], 1, function(x) {
  da <- paste(disorders[grep("X", x)], collapse = ";")
  return(da)
})
data <- tibble::add_column(data, "disorder_association" = anno, .after = "msEntrez")

# Separate rows with multiple disorder associations.
data <- tidyr::separate_rows(data, disorder_association, sep = ";")

# Fix missing disorder annotation for DYRK1A (Ruaud et al., 2015; ID/DD).
idx <- data$Gene == "DYRK1A" & data$disorder_association == ""
data$disorder_association[idx] <- "ID/DD"

# Drop unnecessary columns.
data[, (disorders) := NULL]

# Save the data table.
myfile <- file.path(tabsdir, paste0(script, ".csv"))
fwrite(data, myfile)

# Split into disorder groups.
disorders <- unique(data$disorder_association)
data_list <- data %>%
  group_by(disorder_association) %>%
  group_split()
names(data_list) <- disorders

# Status report.
nGenes <- length(unique(data$msEntrez))
nDisorders <- length(unique(data$disorder_association))
message(paste0(
  "Compiled ", nGenes, " mouse genes associated with ",
  nDisorders, " DBDs!\n"
))

# Check disorder group sizes.
sizes <- sapply(data_list, function(x) length(unique(x$msEntrez)))
message("\nSummary of gene-disease associations:")
knitr::kable(t(sizes), row.names = FALSE, format = "markdown")
message("\n")

# Write as gmt file.
gene_list <- lapply(data_list, function(x) x$msEntrez)
gmt_file <- file.path(gmtdir, paste0(script, ".gmt"))
gmt_source <- "http://dbd.geisingeradmi.org/#additional-information"
write_gmt(gene_list, gmt_source, gmt_file)

# Generate dataset documentation.
documentDataset(gmt_file, short_name,
  Rdir = file.path(root, "R"), datadir
)
