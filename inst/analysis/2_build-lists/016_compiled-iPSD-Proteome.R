#!/usr/bin/env Rscript

## Build a compiled iPSD proteome from several studies:
# [1] Heller et al., 2012
# [2] Kang et al., 2014
# [3] Loh et al., 2016
# [4] Nakamura et al., 2016
# [5] Uezu et al., 2016

## Parameters:
script <- "016_compiled-iPSD-Proteome"
short_name <- "ciPSD"
papers <- list(
  heller = "https://www.ncbi.nlm.nih.gov/pubmed/22768092",
  kang = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4200284/",
  loh = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5167540/",
  nakamura = "https://www.ncbi.nlm.nih.gov/pubmed/27044742",
  uezu = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5432043/"
)

# Load renv
renv::load(getrd())

# Imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(getPPIs)
  library(data.table)
})

# Functions.
suppressWarnings({
  invisible({
    devtools::load_all()
  })
})

# Directories.
root <- getrd()
gmtdir <- file.path(root, "gmt")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "downloads")

# Load the data.
myfile <- list.files(downdir, "iPSD", full.names = TRUE)
data <- read_excel(myfile, sheet = "All Studies")

# Map MGI IDs.
mgi <- data$"MGI:ID"
entrez <- getIDs(mgi, from = "mgi", to = "entrez", species = "mouse")

# Add to data.
data$entrez <- entrez

# Remove na.
data <- subset(data, !is.na(data$entrez))

# Split into groups.
data_list <- data %>%
  group_by(Study) %>%
  group_split()
names(data_list) <- sapply(data_list, function(x) unique(x$Study))

# Coerce to gmt format.
gene_list <- lapply(data_list, function(x) unique(x$entrez))
all_genes <- unique(unlist(gene_list))
names(gene_list) <- gsub(" et al., ", "", names(gene_list))
gene_list[["ciPSD"]] <- all_genes

# Summary.
nGenes <- length(unique(data$entrez))
nStudies <- length(unique(data$Study))
message(paste(
  "Number of iPSD genes compiled from",
  nStudies, "studies:", nGenes
))
df <- sapply(data_list, function(x) length(unique(x$entrez)))
df <- c(df, c("All Studies" = length(all_genes)))
knitr::kable(t(df))

# Save gmt.
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, paste(papers, collapse = ";"), myfile)

# Save as rda and generate documentation.
documentDataset(myfile, short_name, Rdir = file.path(root, "R"), datadir)
