#!/usr/bin/env Rscript

## ARGS:
script <- "056_Gao-IN-PSD95-BioID"
short_name <- "gao2020"
tags <- c("BioID", "proteomics", "interneuron", "excitatory postsynapse")
ref_url <- "in/preparation"
data_url <- "downloads/IN_BioID.xlsx"

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

# Load the data.
message("\nLoading the data from file.")
myfile <- file.path(root, data_url)
all_data <- read_excel_sheets(myfile)

# Collect uniprot ids from the data.
extract_IDs <- function(x) {
  ids <- x %>%
    filter(FDR < 0.1 & logFC > 0) %>%
    dplyr::select(Accession) %>%
    unlist() %>%
    unique()
  return(ids)
}
uniprot_list <- lapply(all_data, extract_IDs)
uniprot <- unique(unlist(uniprot_list))

# Map uniprot to entrez.
message("\nMapping UniprotKB accession to Entrez.")
entrez <- getPPIs::getIDs(uniprot, from = "uniprot", to = "entrez", species = "mouse")
missing <- is.na(entrez)
entrez[missing]

# Map missing ids by hand.
message("\nMapping missing ids by hand using MGI website.")
missing <- is.na(entrez)
# entrez[missing]
mapped_by_hand <- c(
  "P47802" = 17827,
  "Q8VHW2" = 81905,
  "O08585" = 12757,
  "Q61315" = 11789,
  "Q9WTS5" = 23964,
  "Q62083" = 18693,
  "P39087" = 14806,
  "Q8BK63" = 93687,
  "Q8K012" = 214459,
  "Q8CC35" = 104027,
  "P50136" = 12039,
  "Q8BLK3" = 268890,
  "P52189" = 16520,
  "P59823" = 331461,
  "P53026" = 19896,
  "Q8K3E5" = 52906,
  "Q8CIQ7" = 208869
)
entrez[missing] <- mapped_by_hand[names(entrez[missing])]

if (any(is.na(entrez))) {
  stop("Trouble mapping gene names to entrez ids.")
}
if (length(unique(entrez)) != length(entrez)) {
  stop("duplicates!")
}

# Collect list of entrez gene ids.
gene_list <- lapply(uniprot_list, function(x) as.numeric(unique(entrez[x])))

# Status.
ngenes <- length(entrez)
message(paste("\nCompiled", ngenes, "mouse proteins in the PV/CAMKIIA/SST DLG4-BioID proteomes."))
knitr::kable(t(sapply(gene_list, length)))

# Save as gmt, and then save as rda and generate documentation.
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, ref_url, myfile)
documentDataset(myfile, short_name, Rdir, datadir)
