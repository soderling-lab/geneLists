#!/usr/bin/env Rscript

## Scrape LopitDC subcellular marker genes from Geladaki et al., 2019.

## Parameters:
script <- "044_Orre-2019-Subcellular-Proteome"
short_name <- "orre2019"
tags <- c("TMT", "proteomics", "cell line", "subcellular", "proteome")
ref_url <- "https://www.ncbi.nlm.nih.gov/pubmed/30609389"
data_url <- "https://www.cell.com/cms/10.1016/j.molcel.2018.11.035/attachment/9f73533f-b091-4dfa-ae87-ab68c51d9f4a/mmc4.xlsx"

# Load renv.
renv::load(getrd())

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
root <- getrd()
gmtdir <- file.path(root, "gmt")
datadir <- file.path(root, "data")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "downloads")

# Download the data.
message("Download supplemental data from Orre et al., 2019...")
destfile <- file.path(downdir, basename(data_url))
download.file(data_url, destfile, quiet = TRUE)

# Load the data.
all_data <- read_excel_sheets(destfile)

# We will work data from five cell lines:
cell_lines <- c("A431", "MCF7", "H322", "U251", "HCC827")
data <- all_data[cell_lines]

# We will collect protein 'neighborhood' annotations.
# Collect first two columns.
data <- lapply(data, function(x) x[, c(1:2)])

# Stack into a single df.
colNames <- colnames(data[[1]])
df <- dplyr::bind_rows(data, .id = "column_label")
colnames(df)[1] <- "Cell_Line"
colnames(df)[3] <- "Subcellular_Prediction"

# Map Protein gene names to entrez ids.
entrez <- getIDs(df$Protein, from = "symbol", to = "entrez", species = "human")

# How many were not mapped?
not_mapped <- 100 * (sum(is.na(entrez)) / length(entrez))
message(paste("Percent not mapped:", round(not_mapped, 3)))

# Add to data and drop NA.
df$Entrez <- entrez
subdf <- df %>% filter(!is.na(Entrez))

# Map to mouse homologs.
musEntrez <- getHomologs(subdf$Entrez, species = "mouse")
subdf$musEntrez <- musEntrez
subdf <- subdf %>% filter(!is.na(musEntrez))

# Summarize the different fractions:
subdf %>%
  group_by(Subcellular_Prediction) %>%
  summarize(N = length(unique(musEntrez)))

# Build concensus list.
concensus_df <- subdf %>%
  group_by(musEntrez) %>%
  summarize(
    nPred = length(unique(Subcellular_Prediction)),
    Predictions = paste(unique(Subcellular_Prediction), collapse = ";")
  )
concensus_df <- concensus_df %>% filter(nPred == 1)
concensus_list <- split(concensus_df$musEntrez, concensus_df$Predictions)

# Collect list of genes by cell type.
data_list <- subdf %>%
  group_by(Cell_Line) %>%
  group_split()
names(data_list) <- sapply(data_list, function(x) unique(x$Cell_Line))

# For each, collect list of genes.
gene_lists <- lapply(data_list, function(x) {
  split(x$musEntrez, x$Subcellular_Prediction)
})

# Add concensus list.
gene_lists[["concensus"]] <- concensus_list

# For each:
# Save as gmt, and then save as rda and generate documentation.
for (namen in names(gene_lists)) {
  myfile <- file.path(gmtdir, paste(script, namen, sep = "-"))
  write_gmt(gene_lists[[namen]], ref_url, myfile)
  documentDataset(myfile, paste0(short_name, namen),
    Rdir = file.path(root, "R"), datadir
  )
}
