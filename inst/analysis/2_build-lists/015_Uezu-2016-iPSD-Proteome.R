#!/usr/bin/env Rscript

## Scrape iPSD proteome from Uezu et al., 2016
here <- getwd()
root <- dirname(dirname(here))
renv::load(root)

# Urls to data.
data_source <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5432043/"
data_url <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5432043/bin/NIHMS855807-supplement-Table_S4.xlsx"
script <- "015_Uezu-2016-iPSD-Proteome"
short_name <- "iPSD"

# Imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(getPPIs)
  library(data.table)
})

# Functions.
suppressWarnings({
  devtools::load_all()
})

# Directories.
gmtdir <- file.path(root, "gmt")
datadir <- file.path(root, "data")

# Download the data.
myfile <- basename(data_url)
download.file(data_url, myfile, quiet = TRUE)

# Read the data.
data <- read_excel_sheets(myfile, skip = 1)
data <- data[[1]] # There is only one sheet.

# Remove temporary file.
invisible(unlink(myfile))

# Map Uniprot to Entrez.
uniprot <- trimws(data$"UniProt accession")
entrez <- getIDs(uniprot, from = "uniprot", to = "entrez", species = "mouse")

# Use gene symbols to try mapping unmapped genes.
not_mapped <- is.na(entrez)
symbols <- data$"Gene Name"
entrez[not_mapped] <- getIDs(symbols[not_mapped],
  from = "symbol", to = "entrez", species = "mouse"
)

# Map remainder by hand.
not_mapped <- is.na(entrez)
mapped_by_hand <- c("Q8CI71" = 73288, "Q3USH1" = 627214)
entrez[not_mapped] <- mapped_by_hand[names(entrez[not_mapped])]

# One unmapped id is caused by a trailing white space, ignore it.
# sum(is.na(entrez))

# Add to the data.
data <- tibble::add_column(data, entrez, .after = "UniprotID")

# Generate gene groups.
data_list <- data %>%
  group_by(Bait) %>%
  group_split()
names(data_list) <- sapply(data_list, function(x) unique(x$Bait))
gene_groups <- lapply(data_list, function(x) unique(x$entrez))
gene_groups[["iPSD"]] <- unique(data$entrez)

# Status.
nGenes <- length(unique(data$entrez))
message(paste("Collected", nGenes, "iPSD genes."))
sizes <- sapply(gene_groups, length)
knitr::kable(t(sizes), format = "markdown")

# Save as gmt.
ipsd_genes <- gene_groups
myfile <- file.path(gmtdir, script)
write_gmt(ipsd_genes, data_source, myfile)

# Save as rda and generate documentation.
documentDataset(myfile, short_name, Rdir = file.path(root, "R"), datadir)
