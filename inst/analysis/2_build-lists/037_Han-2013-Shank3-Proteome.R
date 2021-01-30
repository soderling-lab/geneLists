#!/usr/bin/env Rscript

## Scrape Han et al., 2013 Shank3 interactome.

## Parameters:
short_name <- "han2013shank3proteome"
script <- "037_Han-2013-Shank3-Proteome"
ref_url <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3923348/"
data_url <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3923348/bin/NIHMS521195-supplement-4.pdf"
pmid <- "24153177"

## Data were converted from pdf to excel and save in /downloads.

renv::load(getrd())

# Imports.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(getPPIs)
  library(readxl)
})

# Functions.
devtools::load_all()

# Directories.
root <- getrd()
gmtdir <- file.path(root, "datasets")
datadir <- file.path(root, "data")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "downloads")

# Load the data.
myfile <- file.path(downdir, "Han-2013-Shank3-Interactome.xlsx")
data <- read_excel(myfile)

# Map mouse symbols to entrez.
entrez <- getIDs(data$"Mouse Symbol",
  from = "symbol", to = "entrez", species = "mouse"
)

# Map missing ids with human symbols.
missing <- is.na(entrez)
hsEntrez <- getIDs(data$"Human Symbol"[missing],
  from = "symbol", to = "entrez", species = "human"
)
entrez[missing] <- getHomologs(hsEntrez, species = "mouse")

# Add to data.
data$entrez <- entrez

# Remove na.
data <- data[!is.na(data$entrez), ]

# Collect lists of genes.
gene_list <- list()
gene_list[["Y2H"]] <- data %>% filter(data$"Y2H" == "o")
gene_list[["IP"]] <- data %>% filter(data$"In vivo IP" == "o")
gene_list <- lapply(gene_list, function(x) unique(x$entrez))

# Summary
df <- data.frame("N Genes" = sapply(gene_list, length))
df <- tibble::add_column(df, "Method" = rownames(df), .before = 1)
knitr::kable(df, row.names = FALSE)

# Save as gene list.
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, ref_url, myfile)

# Save as rda and generate documentation.
documentDataset(myfile, short_name, Rdir = file.path(root, "R"), datadir)
