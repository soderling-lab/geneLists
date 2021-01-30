#!/usr/bin/env Rscript

## Input Parameters:
script <- "049_MS-gene"
short_name <- "msgene"
tags <- c("MS", "multiple sclerosis")
ref_url <- "http://www.msgene.org/"
data_url <- "http://www.msgene.org/"

## CITATION:
# Lill CM, Roehr JT, McQueen MB, Bagade S, Schjeide BM, Zipp F,
# Bertram L. The MSGene Database. Alzheimer Research Forum.
# Available at http://www.msgene.org/. Accessed on 06/14/2020.

# Data were  manually cut-and-paste from the website into excel.
data_file <- "MS-gene.csv" # in root/downloads

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

# Remove row cooresponding to large gene loci.
idx <- grepl("GWA_", data$Gene)
subdat <- subset(data, !idx)

# Map human genes to entrez.
subdat$Entrez <- getIDs(subdat$Gene, from = "Symbol", to = "Entrez", species = "human")

# Map to mouse homologs.
subdat$msEntrez <- getHomologs(subdat$Entrez, species = "mouse")

# Drop un-mapped genes.
mygenes <- unlist(subdat[!is.na(subdat$msEntrez), "msEntrez"]) %>% unique()
gene_list <- list("MS" = mygenes)

# Status.
message(paste(
  "Compiled", length(mygenes), "genes associated with",
  "human Multiple Sclerosis from MSgene.org."
))

# Save as gmt, and then save as rda and generate documentation.
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, ref_url, myfile)
documentDataset(myfile, short_name, Rdir, datadir)
