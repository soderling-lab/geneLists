#!/usr/bin/env Rscript

## Scrape LopitDC subcellular marker genes from Geladaki et al., 2019.

## Input Parameters:
script <- "045_MinoCarta2"
short_name <- "mitocarta2"
tags <- c("mito", "mitochondria", "proteome", "mouse")
ref_url <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4702768/"
data_url <- "http://www.broadinstitute.org/ftp/distribution/metabolic/papers/Pagliarini/MitoCarta2.0/Mouse.MitoCarta2.0.xls"

#--------------------------------------------------------------------

# Load renv.
root <- getrd()
renv::load(root, quiet = TRUE)

# Imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(getPPIs)
  library(TBmiscr)
  library(data.table)
})

# Load functions in root/R.
load_all()

# Directories.
Rdir <- file.path(root, "R")
gmtdir <- file.path(root, "gmt")
datadir <- file.path(root, "data")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "downloads")

# Download the data.
# message("Download supplemental data from Orre et al., 2019...")
destfile <- file.path(downdir, basename(data_url))
download.file(data_url, destfile, quiet = TRUE)

# Load the data.
mysheet <- excel_sheets(destfile)[2] # "A Mouse MitoCarta2.0"
data <- read_xls(destfile, sheet = mysheet)

# All mito carta proteins.
entrez <- data[["MouseGeneID"]]
gene_list <- list("mitocarta2" = entrez)

# Save as gmt, and then save as rda and generate documentation.
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, ref_url, myfile)
documentDataset(myfile, short_name, Rdir, datadir)
