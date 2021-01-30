#!/usr/bin/env Rscript

## ARGS:
script <- "055_Dube-2020-Syp-Proteome"
short_name <- "dube2020"
tags <- c("BioID", "proteomics", "Synapsin", "presynapse")
ref_url <- "in/preparation"
data_url <- "downloads/Synapsin_BioID.xlsx"

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
gmtdir <- file.path(root, "gmt")
datadir <- file.path(root, "data")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "downloads")

# Load the data.
message("\nLoading the data from file.")
myfile <- file.path(root, data_url)
all_data <- read_excel_sheets(myfile)

# Collect uniprot ids and map to entrez.
message("\nMapping UniprotKB accession to Entrez using MGI database.")
syp_data <- all_data[["Final"]]
uniprot <- unique(syp_data$Accession)
entrez <- getPPIs::getIDs(uniprot, from = "uniprot", to = "entrez", species = "mouse")

# Map missing ids by hand.
message("\nMapping missing ids by hand using MGI website.")
missing <- is.na(entrez)
# entrez[missing]
mapped_by_hand <- c(
  "Q4JIM5" = 11352,
  "P47753" = 12340,
  "O54901" = 17470,
  "Q99KN9" = 216705,
  "Q9D8Y0" = 27984,
  "Q9Z0R6" = 20403,
  "Q8BM65" = 241134,
  "O88643" = 18479,
  "Q62083" = 18693,
  "Q9CY18" = 76561,
  "Q8CHC4" = 104015,
  "Q8CC35" = 104027,
  "P0C7L0" = 330319
)

entrez[missing] <- mapped_by_hand[names(entrez[missing])]

if (any(is.na(entrez))) {
  stop("Trouble mapping gene names to entrez ids.")
}
if (length(unique(entrez)) != length(entrez)) {
  stop("duplicates!")
}

# Status.
ngenes <- length(entrez)
message(paste("\nCompiled", ngenes, "mouse proteins in the SYP-iBioID presynaptic proteome."))

# Save this as a gene list.
gene_list <- list("syp-presynapse" = entrez)

# Save as gmt, and then save as rda and generate documentation.
Rdir <- file.path(root, "R")
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, ref_url, myfile)
documentDataset(myfile, short_name, Rdir, datadir)
