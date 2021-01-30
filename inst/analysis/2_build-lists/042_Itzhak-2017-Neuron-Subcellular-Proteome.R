#!/usr/bin/env Rscript

## Scrape LopitDC subcellular marker genes from Geladaki et al., 2019.

## Parameters:
script <- "042_Itzhak-2017-Neuron_Subcellular-Proteome"
short_name <- "itzhak2017"
tags <- c("TMT", "mouse primary neurons", "neurons", "subcellular", "proteome")
ref_url <- "https://www.ncbi.nlm.nih.gov/pubmed/28903049"
data_urls <- c(
  S1 = "https://bit.ly/3b5cdg5",
  S3 = "https://bit.ly/3cio2Qv",
  S4 = "https://bit.ly/3abXE9z"
)

renv::load(getrd())

# Imports.
suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
  library(dplyr)
  library(getPPIs)
})

# Functions.
suppressWarnings({
  devtools::load_all()
})

# Directories.
root <- getrd()
gmtdir <- file.path(root, "gmt")
datadir <- file.path(root, "data")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "downloads")

# Download the data.
message("Download supplemental files from Itzhak et al., 2017...")
destfiles <- file.path(downdir, paste0(names(data_urls), ".xlsx"))
names(destfiles) <- names(data_urls)
invisible({
  mapply(download.file, data_urls, destfiles, quiet = TRUE)
})

# Load the data from S4.
all_data <- read_excel_sheets(destfiles["S4"])
data <- all_data[[1]]

# Clean-up column names.
colnames(data) <- gsub("\\ ", "_", colnames(data))

# Drop rows in which there was no compartment prediction.
data <- data %>% filter(Prediction != "NA")

# Get Uniprot IDs.
uniprot <- data$Canonical_lead_protein_ID

# Map mouse uniprot IDs to entrez.
entrez <- getIDs(uniprot, from = "uniprot", to = "entrez", species = "mouse")

# Map missing entrez using gene symbols.
not_mapped <- is.na(entrez)
symbol <- data$Lead_gene_name[not_mapped]
entrez[not_mapped] <- getIDs(symbol, from = "symbol", to = "entrez", species = "mouse")

# Check any remaining unmapped ids.
not_mapped <- is.na(entrez)
entrez[not_mapped]

# Map remaining missing by hand.
mapped_by_hand <- c(
  99633, 110329271, 243780, NA, 18563, 14467, 68195,
  109135, 67532, 60365, NA, NA, 338320, NA,
  NA
)
entrez[not_mapped] <- mapped_by_hand
data <- tibble::add_column(data,
  "Entrez" = entrez,
  .after = "Canonical_lead_protein_ID"
)

# Check how many are not mapped to mouse homologs.
not_mapped <- sum(is.na(data$Entrez))

# Remove unmapped rows.
data <- data %>% filter(Entrez != "NA")

# Collect list of protein subcellular annotations.
pred_list <- split(data$Entrez, data$Prediction)
marker_list <- split(data$Entrez, data$"Marker_protein?")

# Summary:
n <- sapply(pred_list, length)
df <- as.data.frame(n)
df <- tibble::add_column(df,
  "Predicted Subcellular Compartment" = rownames(df),
  .before = 1
)
knitr::kable(df, row.names = FALSE)

# Summary:
n <- sapply(marker_list, length)
df <- as.data.frame(n)
df <- tibble::add_column(df,
  "Subcellular Compartment" = rownames(df),
  .before = 1
)
knitr::kable(df, row.names = FALSE)

# Save as gmt, and then save as rda and generate documentation.
# Predictions:
myfile <- file.path(gmtdir, paste(script, "predictions", sep = "-"))
write_gmt(pred_list, ref_url, myfile)
documentDataset(myfile, paste0(short_name, "predictions"), Rdir = file.path(root, "R"), datadir)

# Markers:
myfile <- file.path(gmtdir, paste(script, "markers", sep = "-"))
write_gmt(marker_list, ref_url, myfile)
documentDataset(myfile, paste0(short_name, "markers"), Rdir = file.path(root, "R"), datadir)
