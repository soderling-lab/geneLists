#!/usr/bin/env Rscript

## Scrape LopitDC subcellular marker genes from Geladaki et al., 2019.

renv::load(getrd())

## Parameters:
script <- "010_Geladaki-2019-HeLa-LopitDC-Proteome"
short_name <- "lopitDC"
tags <- c("TMT", "Lopit", "HyperLopit", "LopitDC", "HyperLopitDC")
ref_url <- "https://www.ncbi.nlm.nih.gov/pubmed/30659192"
urls <- c(
  LopitDC = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018-08191-w/MediaObjects/41467_2018_8191_MOESM4_ESM.xlsx",
  Lopit = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018-08191-w/MediaObjects/41467_2018_8191_MOESM5_ESM.xlsx"
)

# Imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(getPPIs)
  library(data.table)
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
download.files(urls, quiet = TRUE)

# Load the data.
myfiles <- basename(urls)
data <- lapply(myfiles, function(file) {
  data <- lapply(excel_sheets(file), function(sheet) {
    df <- read_excel(file, sheet, skip = 1)
    colnames(df)[1] <- "Uniprot"
    return(df)
  })
  names(data) <- excel_sheets(file)
  return(data)
})
names(data) <- names(urls)

# Remove downloaded files.
unlink(myfiles)

# Loop to collect data of interest
cols <- c("Uniprot", "Subcellular markers", "Curated SVM predictions", "SVM score")
markers <- list()
for (dataset in seq_along(data)) {
  subdat <- data[[dataset]]
  subdat <- subdat[[1]] %>% dplyr::select(all_of(cols))
  markers[[dataset]] <- subdat
}
names(markers) <- names(data)

# Let's focus on the LopitDC data.
markers <- markers[["LopitDC"]]

# Map Uniprot IDs to gene names.
# Insure that we have organism mapping databases installed.
uniprot <- sapply(strsplit(markers$Uniprot, "-"), "[", 1)
entrez <- getIDs(uniprot, from = "uniprot", to = "entrez", species = "human")
markers <- tibble::add_column(markers, "Entrez" = entrez, .after = "Uniprot")

# Remove NA.
markers <- markers[!is.na(markers$Entrez), ]

# Map Human entrez to mouse homologs.
msEntrez <- getHomologs(markers$Entrez, taxid = 10090)
markers <- tibble::add_column(markers, "msEntrez" = msEntrez, .after = "Entrez")

# Remove NA.
markers <- markers[!is.na(markers$msEntrez), ]

# Add mouse gene symbols.
symbols <- getIDs(markers$msEntrez,
  from = "entrez", to = "symbol", species = "mouse"
)
markers <- tibble::add_column(markers, "msSymbol" = symbols, .after = "msEntrez")

# Collect list of genes and their subcellular localization or
# predicted subcellular localization.
idy <- c("Subcellular markers", "Curated SVM predictions")
fractions <- lapply(idy, function(y) split(markers, markers[[y]]))
names(fractions) <- c("Markers", "Predictions")

# Status
sizes <- lapply(fractions, function(x) sapply(x, nrow))

# Markers.
df <- data.frame("Proteins" = sizes[[1]])
df <- tibble::add_column(df, Fraction = rownames(df), .before = 1)
knitr::kable(df, row.names = FALSE)

# Predictions.
df <- data.frame("Proteins" = sizes[[2]])
df <- tibble::add_column(df, Fraction = rownames(df), .before = 1)
knitr::kable(df, row.names = FALSE)

# Save as gene lists.
gene_list1 <- sapply(fractions[[1]], function(x) unique(x$msEntrez))
gene_list2 <- sapply(fractions[[2]], function(x) unique(x$msEntrez))
myfiles <- file.path(gmtdir, paste(script, names(fractions), sep = "_"))
write_gmt(gene_list1, ref_url, myfiles[1]) # Markers
write_gmt(gene_list2, ref_url, myfiles[2]) # Predictions

# Save as rda and generate documentation.
documentDataset(myfiles[1], paste0(short_name, "markers"),
  Rdir = file.path(root, "R"), datadir
)
documentDataset(myfiles[2], paste0(short_name, "predictions"),
  Rdir = file.path(root, "R"), datadir
)
