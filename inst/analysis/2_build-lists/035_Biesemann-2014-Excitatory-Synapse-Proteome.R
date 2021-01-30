#!/usr/bin/env Rscript

## Scrape Synaptic vesicle proteome from Takamori et al., 2006

## Parameters:
short_name <- "biesemann2014sortedsynaptosomes"
script <- "035_Biesemann-2014-Excitatory-Synapse-Proteome"
ref_url <- "http://www.ncbi.nlm.nih.gov/pubmed?term=24413018"
data_url <- c(
  S1 = "https://www.embopress.org/action/downloadSupplement?doi=10.1002%2Fembj.201386120&file=embj201386120-sup-0005-TableS1.xls",
  S2 = "https://www.embopress.org/action/downloadSupplement?doi=10.1002%2Fembj.201386120&file=embj201386120-sup-0006-TableS2.xls",
  S3 = "https://www.embopress.org/action/downloadSupplement?doi=10.1002%2Fembj.201386120&file=embj201386120-sup-0007-TableS3.xls",
  S4 = "https://www.embopress.org/action/downloadSupplement?doi=10.1002%2Fembj.201386120&file=embj201386120-sup-0008-TableS4.xlsx"
)

# Load renv:
renv::load(getrd())

# Imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(getPPIs)
  library(data.table)
})

# Functions.
devtools::load_all()

# Directories.
root <- getrd()
gmtdir <- file.path(root, "datasets")
datadir <- file.path(root, "data")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "downloads")

# Download the data.
myfiles <- file.path(downdir, basename(data_url))
download.files(data_url, destdir = downdir, quiet = TRUE)

## Read the data.
# S1 is Synaptosome.
# S2 is FASS - sorted synaptosomes
# S3 is scaffold file.
# S4 is enriched/depleted proteins.
all_data <- lapply(myfiles, read_excel)
names(all_data) <- paste0("S", c(1:4))

# Loop to map gi ids using uniprot API with python script.
for (i in seq(length(all_data))) {
  data <- all_data[[i]]
  idy <- grep("gi-", colnames(data))
  ids <- as.character(data[[idy]])
  if (any(grepl("\\|", ids))) {
    ids <- sapply(strsplit(ids, "\\|"), "[", 2)
  }
  IDs <- file.path(getwd(), "IDs.txt")
  writeLines(ids, IDs)
  # Call python script to map identifiers
  pyScript <- file.path(root, "Py", "mapIds.py")
  cmd <- paste(pyScript, IDs, "P_GI", "P_ENTREZGENEID")
  result <- system(cmd, intern = TRUE)
  # Remove temporary file.
  unlink(IDs)
  # Parse the result.
  entrez <- strsplit(gsub("\\[|\\]|'", "", result), ", ")[[1]]
  entrez[entrez == "None"] <- NA
  # Status.
  not_mapped <- 100 * sum(is.na(entrez)) / length(entrez)
  message(paste(
    "Percent successfully mapped genes:",
    round(100 - not_mapped, 3), "%."
  ))
  # Map to mouse.
  msEntrez <- getHomologs(entrez, species = "mouse")
  data <- tibble::add_column(data, entrez, .after = idy)
  data <- tibble::add_column(data, msEntrez, .after = idy + 1)
  # Drop NA, write over original df.
  data <- data[!is.na(data$msEntrez), ]
  all_data[[i]] <- data
}

## Create a list of genes.
gene_list <- list()
gene_list[["All Synaptosome"]] <- unique(all_data[[1]]$msEntrez)
gene_list[["All FASS"]] <- unique(all_data[[2]]$msEntrez)
gene_list[["FASS Enriched (>2)"]] <- unique(all_data[[4]]$msEntrez)

# Summary
df <- data.frame("Proteins" = sapply(gene_list, length))
df <- tibble::add_column(df, "Protein Class" = rownames(df), .before = 1)
knitr::kable(df, row.names = FALSE)

# Save as gene list.
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, ref_url, myfile)

# Save as rda and generate documentation.
documentDataset(myfile, short_name, Rdir = file.path(root, "R"), datadir)
