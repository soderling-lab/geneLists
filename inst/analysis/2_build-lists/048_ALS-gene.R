#!/usr/bin/env Rscript

## Input Parameters:
script <- "048_ALS-gene"
short_name <- "alsgene"
tags <- c("ALS")
ref_url <- "https://ghr.nlm.nih.gov/gene"
data_url <- "https://ghr.nlm.nih.gov/gene"

# Data were  manually cut-and-paste from the website into excel.
data_file <- "ALS-gene.csv" # in root/downloads

# Load renv.
root <- getrd()
renv::load(root, quiet = TRUE)

# Imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(getPPIs)
  library(data.table)
})

# Load functions in root/R.
suppressMessages({
  devtools::load_all()
})

# Directories.
Rdir <- file.path(root, "R")
gmtdir <- file.path(root, "gmt")
datadir <- file.path(root, "data")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "downloads")

# Load the data.
myfile <- file.path(downdir, data_file)
data <- fread(myfile)

# Map human genes to entrez.
data$Entrez <- getIDs(data$Gene, from = "Symbol", to = "Entrez", species = "human")

# Map to mouse homologs.
data$msEntrez <- getHomologs(data$Entrez, species = "mouse")

# Drop un-mapped genes.
mygenes <- unlist(data[!is.na(data$msEntrez), "msEntrez"]) %>% unique()
gene_list <- list("ALS" = mygenes)

# Status.
message(paste(
  "Compiled", length(mygenes), "genes associated with",
  "human ALS from ALSgene.org"
))

# Save as gmt, and then save as rda and generate documentation.
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, ref_url, myfile)
documentDataset(myfile, short_name, Rdir, datadir)
