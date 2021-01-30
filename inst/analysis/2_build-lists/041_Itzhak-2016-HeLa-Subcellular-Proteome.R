#!/usr/bin/env Rscript

## Scrape subcellular marker protein from Itzhak et al., 2016

## Parameters:
script <- "041_Itzhak-2016-HeLa-Subcellular-Proteome"
short_name <- "itzhak2016"
tags <- c("HeLa", "subcellular fractionation", "proteomics")
ref_url <- "https://www.ncbi.nlm.nih.gov/pubmed/27278775"

# Urls to the data, shortened with twzer, see bin/shorten.
data_urls <- c(
  S1 = "https://bit.ly/2V54pFG",
  S4 = "https://bit.ly/3b8GZ7O",
  S5 = "https://bit.ly/3a6QpQ1"
)

renv::load(getrd())

# Imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(getPPIs)
  library(data.table)
})

# Load functions.
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
destfiles <- file.path(downdir, paste0(names(data_urls), ".xlsx"))
names(destfiles) <- names(data_urls)
invisible({
  mapply(download.file, data_urls, destfiles, quiet = TRUE)
})

# Load the data in S1 -- HeLa spatial protoeme.
sheets <- excel_sheets(destfiles["S1"])
all_data <- read_excel_sheets(destfiles["S1"])
data <- all_data$"Compact HeLa Spatial Proteome"

# Clean-up column names.
colnames(data) <- gsub("-", "", colnames(data))
colnames(data) <- gsub("\\ ", "_", colnames(data))

# Drop rows in which there was no compartment prediction.
data <- data %>% filter(Compartment_Prediction != "No Prediction")

# Get Uniprot IDs.
uniprot <- data$"Lead_Protein_ID"

# Remove isoform specifier.
uniprot <- gsub("-[1-9]{1,2}", "", uniprot)

# Map human uniprot IDs to entrez.
entrez <- getIDs(uniprot, from = "uniprot", to = "entrez", species = "human")

# Map missing entrez using gene symbols.
not_mapped <- is.na(entrez)
symbol <- data$Lead_Gene_name[not_mapped]
entrez[not_mapped] <- getIDs(symbol, from = "symbol", to = "entrez", species = "human")

# Check any remaining unmapped ids.
not_mapped <- is.na(entrez)
# entrez[not_mapped]

# Map by hand.
entrez[not_mapped] <- c(10207, 112268437, 23392)
data <- tibble::add_column(data, "Entrez" = entrez, .after = "Lead_Protein_ID")

# Map Human entrez to mouse homologs.
msEntrez <- getHomologs(data$Entrez, taxid = 10090)
data <- tibble::add_column(data, "msEntrez" = msEntrez, .after = "Entrez")

# Check how many are not mapped to mouse homologs.
not_mapped <- sum(is.na(data$msEntrez))

# Remove unmapped rows.
data <- data %>% filter(msEntrez != "NA")

# Collect list of annotations.
pred_list <- split(data$msEntrez, data$Compartment_Prediction)

# Summary:
message("Summary of the composition of prediction subcellular compartments:")
n <- sapply(pred_list, length)
df <- as.data.frame(n)
df <- tibble::add_column(df, "Subcellular Compartment" = rownames(df), .before = 1)
knitr::kable(df, row.names = FALSE)

# Save as gene list.
myfile <- file.path(gmtdir, paste(script, "predictions", sep = "-"))
write_gmt(pred_list, ref_url, myfile) # Markers

# Save as rda and generate documentation.
documentDataset(myfile, paste0(short_name, "predictions"), Rdir = file.path(root, "R"), datadir)

## Repeat this process for the HeLa cell markers:
# Work on the markers.
markers <- all_data[["Organellar Markers HeLa"]]
colnames(markers) <- gsub("-", "", colnames(markers))
colnames(markers) <- gsub("\\ ", "_", colnames(markers))
uniprot <- markers$"Protein_ID_(canonical)"
entrez <- getIDs(uniprot, from = "uniprot", to = "entrez", species = "human")
not_mapped <- is.na(entrez)
entrez["Q9BRJ2"] <- 84311
markers <- tibble::add_column(markers, Entrez = entrez, .after = "Gene_name")
msEntrez <- getHomologs(markers$Entrez, taxid = 10090)
markers <- tibble::add_column(markers, "msEntrez" = msEntrez, .after = "Entrez")
markers <- markers %>% filter(msEntrez != "NA")
marker_list <- split(markers$msEntrez, markers$Compartment)

# Summary:
message("Summary of subcellular compartment markers:")
n <- sapply(marker_list, length)
df <- as.data.frame(n)
df <- tibble::add_column(df, "Subcellular Compartment" = rownames(df), .before = 1)
knitr::kable(df, row.names = FALSE)

# Save as gene list.
myfile <- file.path(gmtdir, paste(script, "markers", sep = "-"))
write_gmt(marker_list, ref_url, myfile)

# Save as rda and generate documentation.
documentDataset(myfile, paste0(short_name, "markers"), Rdir = file.path(root, "R"), datadir)
