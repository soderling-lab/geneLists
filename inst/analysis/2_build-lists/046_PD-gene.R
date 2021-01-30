#!/usr/bin/env Rscript

## Input Parameters:
script <- "046_PD-gene"
short_name <- "pdgene"
tags <- c("PD", "parkinsons disease", "genes")
ref_url <- "https://pubmed.ncbi.nlm.nih.gov/25064009/"
data_url <- "http://www.pdgene.org/top_results"

# Data were  manually cut-and-paste from the website into excel.
data_file <- "PD-gene.csv" # in root/downloads

#--------------------------------------------------------------------

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

# Remove intergenic rows.
subdat <- data %>% filter(Gene != "intergenic")

# Map human genes to entrez.
Entrez <- getIDs(subdat$Gene, from = "Symbol", to = "Entrez", species = "human")
subdat <- tibble::add_column(subdat, Entrez, .after = "Gene")

# Map to mouse homologs.
msEntrez <- getHomologs(subdat$Entrez, species = "mouse")
subdat <- tibble::add_column(subdat, msEntrez, .after = "Entrez")

# Drop un-mapped genes.
mygenes <- subdat %>%
  filter(msEntrez != "NA") %>%
  dplyr::select(msEntrez) %>%
  unlist() %>%
  unique()
gene_list <- list("PD" = mygenes)

# Status.
message(paste(
  "Compiled", length(mygenes), "genes associated with",
  "human Parkinson's disease."
))

# Save as gmt, and then save as rda and generate documentation.
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, ref_url, myfile)
documentDataset(myfile, short_name, Rdir, datadir)
