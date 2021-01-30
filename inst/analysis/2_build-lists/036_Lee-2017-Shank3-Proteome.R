#!/usr/bin/env Rscript

## Scrape Lee et al., 2017 Shank3 interactome.

## Parameters:
short_name <- "lee2017shank3proteome"
script <- "036_Lee-2017-Shank3-Proteome"
ref_url <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5395616/"
data_url <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5395616/bin/Table_1.XLSX"

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

# Download the data.
myfile <- file.path(downdir, basename(data_url))
download.file(data_url, myfile, quiet = TRUE)

# Load the data into R.
data <- read_excel(myfile, sheet = 1, skip = 1)

# Map mgi to entrez.
data$entrez <- getIDs(data$MGI, from = "mgi", to = "entrez", species = "mouse")

# Collect lists of genes.
gene_list <- list()
gene_list[["mPFC"]] <- data %>% filter(data$"mPFC" == "Y")
gene_list[["HP+STR"]] <- data %>% filter(data$"HP+STR" == "Y")
gene_list[["Y2H"]] <- data %>% filter(data$"Y2H" == "Y")
gene_list <- lapply(gene_list, function(x) unique(x$entrez))

# Summary
df <- data.frame("N Genes" = sapply(gene_list, length))
df <- tibble::add_column(df, "Class" = rownames(df), .before = 1)
knitr::kable(df, row.names = FALSE)

# Save as gene list.
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, ref_url, myfile)

# Save as rda and generate documentation.
documentDataset(myfile, short_name, Rdir = file.path(root, "R"), datadir)
