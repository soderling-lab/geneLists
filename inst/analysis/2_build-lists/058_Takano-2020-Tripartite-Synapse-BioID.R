#!/usr/bin/env Rscript

## ARGS:
script <- "058_Takano-2020-Tripartite-Synapse-BioID"
short_name <- "takano2020"
tags <- c("BioID", "proteomics", "synapse", "astrocyte", "tripartite")
ref_url <- "in/preparation"
data_url <- "downloads/Tripartite_BioID.csv"

# Load renv.
here <- getwd()
root <- dirname(dirname(here))
renv::load(root)

# Imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(getPPIs)
  library(TBmiscr)
  library(data.table)
})

# Load functions in root/R.
TBmiscr::load_all()

# Directories.
root <- getrd()
Rdir <- file.path(root, "R")
gmtdir <- file.path(root, "gmt")
datadir <- file.path(root, "data")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "downloads")

# Load the data from root/downloads.
message("\nLoading the data from file.")
myfile <- file.path(root, data_url)
df <- fread(myfile)

# Extract entrez ids.
entrez <- df$"Entrez Gene ID"

# Save as gene list.
gene_list <- list('takano2020' = entrez)

# Save as gmt, and then save as rda and generate documentation.
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, ref_url, myfile)
documentDataset(myfile, short_name, Rdir, datadir)
