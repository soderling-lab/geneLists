#!/usr/bin/env Rscript

## Scrape LopitDC subcellular marker genes from Geladaki et al., 2019.

## Input Parameters:
script <- "057_GSEA-Hallmark-Genes"
short_name <- "hallmark"
tags <- c("gsea", "hallmark", "mouse")
ref_url <- "NA"
data_url <- "/downloads/h.all.v7.1.entrez.gmt"

#--------------------------------------------------------------------
## Prepare the workspace
#--------------------------------------------------------------------

# Load renv.
here <- getwd()
root <- dirname(dirname(here))
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
TBmiscr::load_all()

# Directories.
Rdir <- file.path(root, "R")
gmtdir <- file.path(root, "gmt")
datadir <- file.path(root, "data")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "downloads")

# Load the data.
gene_list <- read_gmt(file.path(root, data_url))

# Map to mouse.
mouse_genes <- list()
message("\nMapping human Hallmark genes to mouse  homologs.")
pbar <- txtProgressBar(max = length(gene_list), style = 3)
for (i in c(1:length(gene_list))) {
  mouse_genes[[i]] <- getPPIs::getHomologs(gene_list[[i]], species = "mouse")
  setTxtProgressBar(pbar, i)
}

# Remove unmapped genes.
names(mouse_genes) <- names(gene_list)
filt_genes <- lapply(mouse_genes, function(x) x[!is.na(x)])

# Status.
message("\nSummary of mouse Hallmark genes:")
knitr::kable(sapply(filt_genes, length))

# Save human data as gmt, and then save as rda and generate documentation.
myfile <- file.path(gmtdir, paste0("hs", script))
write_gmt(gene_list, ref_url, myfile)
documentDataset(myfile, paste0("hs", short_name), Rdir, datadir)

# Save mouse data as gmt, and then save as rda and generate documentation.
myfile <- file.path(gmtdir, paste0("ms", script))
write_gmt(mouse_genes, ref_url, myfile)
documentDataset(myfile, paste0("ms", short_name), Rdir, datadir)

message("\nDone!")
