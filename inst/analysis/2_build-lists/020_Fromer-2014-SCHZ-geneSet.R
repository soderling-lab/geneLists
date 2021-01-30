#!/usr/bin/env Rscript

## Scrape Schizophrenia genes from Fromer et al., 2014

## Parameters:
short_name <- "fromer2014schz"
script <- "020_Fromer-2014-SCHZ-geneSet"
ref_url <- "https://www.ncbi.nlm.nih.gov/pubmed/24463507"
data_url <- c(
  S1 = "https://static-content.springer.com/esm/art%3A10.1038%2Fnature12929/MediaObjects/41586_2014_BFnature12929_MOESM39_ESM.xlsx",
  S2 = "https://static-content.springer.com/esm/art%3A10.1038%2Fnature12929/MediaObjects/41586_2014_BFnature12929_MOESM40_ESM.xlsx"
)

# Load renv:
renv::load(getrd())

# Imports.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(getPPIs)
})

# Functions.
devtools::load_all()

# Directories.
root <- getrd()
gmtdir <- file.path(root, "gmt")
datadir <- file.path(root, "data")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "downloads")

# Download the data.
myfiles <- basename(data_url)
download.files(data_url, quiet = TRUE)

# Load the data.
data <- list()
data[["S1"]] <- read_excel_sheets(myfiles[1])
data[["S2"]] <- read_excel_sheets(myfiles[2])
data <- unlist(data, recursive = FALSE)

# Remove downloaded files.
unlink(myfiles)

# Drop individuals data from list.
n_individuals <- data[[2]]
data[[2]] <- NULL

# Loop to clean up Genes and Gene annotations columns.
data <- lapply(data, function(df) {
  df <- tidyr::separate_rows(df, "Gene annotations", sep = ",")
  df <- tidyr::separate_rows(df, "Genes", sep = ",")
  return(df)
})

# Map human genes to entrez.
data <- lapply(data, function(df) {
  entrez <- getIDs(df$Genes, from = "symbol", to = "entrez", species = "human")
  df <- tibble::add_column(df, entrez, .after = "Genes")
  return(df)
})

# Map human entrez to mouse.
data <- lapply(data, function(df) {
  msEntrez <- getHomologs(df$entrez, species = "mouse")
  df <- tibble::add_column(df, msEntrez, .after = "entrez")
  return(df)
})

# Remove NAs.
data <- lapply(data, na.omit)

# Save table S1.
# myfile <- file.path(tabsdir, paste0(script, ".csv"))
# fwrite(data[[1]], myfile)

# Save table S2.
# myfile <- file.path(tabsdir, paste0(
#  gsub("SCHZ-geneSet", "", script),
#  "Compiled-Literature-DBD-geneSet.csv"
# ))
# fwrite(data[[2]], myfile)

# Summary of mutation types.
df <- data[[1]]
table(df$"Gene annotations")

# Remove silent mutations.
df <- df %>% filter(df$"Gene annotations" != "silent")

# Status
genes <- unique(df$msEntrez)
nGenes <- length(genes)
message(paste(
  "Compiled", nGenes, "mouse genes with de novo",
  "mutations in humans with schizophrenia."
))

# Save as gene list.
myfile <- file.path(gmtdir, script)
gene_list <- list("SCHZ" = genes)
write_gmt(gene_list, ref_url, myfile)

# Save as rda and generate documentation.
documentDataset(myfile, short_name, Rdir = file.path(root, "R"), datadir)

#--------------------------------------------------------------------
## Save the other datasets compiled by Fromer et al.
#--------------------------------------------------------------------

# Remove controls and unaffected individuals.
df <- data[[2]]
out <- c("CONTROL", "Unaffected sibling")
idx <- df$"Child phenotype" %in% out
df <- df[!idx, ]

# Check mutation types.
table(df$"Gene annotations")

# Remove silent mutations.
df <- df %>% filter(df$"Gene annotations" != "silent")

# Remove mutations of unknown consequence.
df <- df %>% filter(df$"Gene annotations" != "exonic-unknown")

# Group data by study.
data_list <- split(df, df$Study)

# Rename data in same format used by other gmt lists.
all_studies <- unique(gsub("\\ ", "_", n_individuals$STUDY))
idx <- sapply(names(data_list), function(x) grep(x, all_studies))
studies <- all_studies[idx]
studies <- gsub("_", "", studies)
studies <- tolower(gsub(",", "", studies))
studies <- paste0(studies, sapply(data_list, function(x) unique(x$"Child phenotype")))
names(data_list) <- studies

# Get urls to studies.
data_urls <- c(
  "https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000298.v4.p3",
  "https://www.ncbi.nlm.nih.gov/pubmed/23033978",
  "https://www.ncbi.nlm.nih.gov/pubmed/22456906",
  "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3894107/",
  "https://www.ncbi.nlm.nih.gov/pubmed/22542183",
  "https://www.ncbi.nlm.nih.gov/pubmed/22495311",
  "https://www.ncbi.nlm.nih.gov/pubmed/23160955",
  "https://www.ncbi.nlm.nih.gov/pubmed/23020937",
  "https://www.ncbi.nlm.nih.gov/pubmed/22495306",
  "https://www.ncbi.nlm.nih.gov/pubmed/23042115"
)
# ARRA = Autism Sequencing Collaboration

# Collect lists of genes.
disease_anno <- c("ASD", "ID", "SZ", "SZ", "ASD", "ASD", "ASD", "ID", "ASD", "SZ")
gene_list <- lapply(data_list, function(x) unique(x$msEntrez))
names(gene_list) <- paste(disease_anno, names(gene_list), sep = ":")

# Summary.
sizes <- sapply(gene_list, length)
tab <- data.table("N Genes" = sapply(gene_list, length))
tab <- tibble::add_column(tab, "Study" = names(sizes), .before = 1)
knitr::kable(tab, row.names = FALSE)

# Save gmt file.
myfile <- file.path(gmtdir, gsub("SCHZ", "DBD", script))
write_gmt(gene_list, paste(data_urls, collapse = ";"), myfile)

# Save as rda and generate documentation.
documentDataset(myfile, "fromer2014dbd", Rdir = file.path(root, "R"), datadir)
