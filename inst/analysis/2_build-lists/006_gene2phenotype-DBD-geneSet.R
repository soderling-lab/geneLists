#!/usr/bin/env Rscript

## DBD genes compiled from G2P database.

renv::load(getrd())

## User parameters:
script <- "006_gene2phenotype-DBD-geneSet"
short_name <- "g2pDBD"
dataset <- "DD" # One of c("Cancer","DD","Eye","Skin")
diseases <- c(
  "autism",
  "intellectual disability",
  "attention deficit hyperactivity disorder",
  "schizophrenia",
  "bipolar disorder",
  "epilepsy"
)

# Imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(getPPIs)
  library(data.table)
})

# Functions.
devtools::load_all()

# Directories.
here <- getwd()
root <- getrd()
gmtdir <- file.path(root, "gmt")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "downloads")

# Downoad the raw data.
base_url <- "https://www.ebi.ac.uk/gene2phenotype/downloads"
datasets <- c(
  Cancer = "CancerG2P",
  DD = "DDG2P",
  Eye = "EyeG2P",
  Skin = "SkinG2P"
)
myurl <- file.path(base_url, paste0(datasets[dataset], ".csv.gz"))
download.file(myurl, destfile = basename(myurl), quiet = TRUE)
system(command = paste("gunzip", basename(myurl)))
myfile <- paste0(datasets[dataset], ".csv")
raw_data <- fread(myfile)
unlink(myfile)

# Fix column names.
colnames(raw_data) <- gsub(" ", "_", colnames(raw_data))

# Map human gene symbols to Entrez.
genes <- unique(raw_data$gene_symbol)
nHsGenes <- length(unique(genes))
entrez <- getIDs(genes, from = "symbol", to = "entrez", species = "human")
names(entrez) <- genes

# Map missing Entrez IDs by hand.
not_mapped <- genes[is.na(entrez)]
mapped_by_manual_search <- c(
  "KIF1BP" = 26128,
  "KARS" = 3735,
  "MT-TP" = 4571,
  "QARS" = 5859,
  "HARS" = 3035,
  "HIST3H3" = 8290,
  "AARS" = 16,
  "ARSE" = 415,
  "DARS" = 1615,
  "HIST1H4B" = 8366,
  "IARS" = 3376,
  "HIST1H4C" = 8364,
  "HIST1H1E" = 3006,
  "HIST1H4J" = 8363,
  "EPRS" = 2058,
  "CARS" = 833,
  "TARS" = 689,
  "HIST1H2AC" = 8334,
  "ADPRS" = 54936
)
entrez[not_mapped] <- mapped_by_manual_search[names(entrez[not_mapped])]

# Check.
check <- sum(is.na(entrez)) == 0
if (!check) {
  stop("Not all genes were mapped to entrez!")
}

# Add entrez IDs to data.
idy <- "gene_symbol" # Column after which Entrez ids will be added.
data <- tibble::add_column(raw_data, "Entrez" = entrez[raw_data[[idy]]], .after = idy)

# Map human genes in to their mouse homologs.
hsEntrez <- data$Entrez
msEntrez <- getHomologs(hsEntrez, taxid = 10090)
data <- tibble::add_column(data, msEntrez = msEntrez, .after = "Entrez")
# Remove rows with unmapped genes.
data <- data[msEntrez != "NA"]

# Get disorders that affect the brain/cognition.
disease_types <- "Brain/Cognition"
data <- tidyr::separate_rows(data, organ_specificity_list, sep = ";")
data <- subset(data, data$organ_specificity_list == disease_types)

# Collect diseases of interest.
rows <- lapply(diseases, function(x) grep(x, tolower(data$disease_name)))

# Save table.
idx <- unlist(rows)
myfile <- file.path(tabsdir, paste0(script, ".csv"))
fwrite(data[idx, ], myfile)

# Split into groups.
data_list <- lapply(rows, function(idx) data[idx, ])
names(data_list) <- diseases

# Check disorder group sizes.
sizes <- sapply(data_list, function(x) length(unique(x$msEntrez)))
message("\nSummary of gene-disease associations:")
knitr::kable(t(sizes), row.names = FALSE, format = "markdown")
message("\n")

# Remove groups with 0 genes.
data_list <- data_list[-which(sizes == 0)]

# Status report.
nGenes <- length(unique(unlist(sapply(data_list, function(x) x$msEntrez))))
nDisorders <- length(data_list)
message(paste("Compiled", nGenes, "mouse genes associated with", nDisorders, "DBDs!"))

# save as gmt.
gmt_list <- lapply(data_list, function(x) x$msEntrez)
gmt_source <- "https://www.ebi.ac.uk/gene2phenotype/downloads"
gmt_file <- file.path(gmtdir, script)
write_gmt(gmt_list, gmt_source, gmt_file)

# Save as rda and generate documentation.
documentDataset(gmt_file, short_name, Rdir = file.path(root, "R"), datadir)
