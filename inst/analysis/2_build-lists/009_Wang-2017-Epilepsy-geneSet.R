#!/usr/bin/env Rscript

## Building an Epilepsy gene dataset from the genes compiled by
# Wang et al., 2017. Data were obtained from the personal communication
# with the author.

renv::load(getrd())

## User parameters:
script <- "009_Wang-2017-Epilepsy-geneSet"
short_name <- "wang2017Epilepsy"
paper_url <- "https://www.ncbi.nlm.nih.gov/pubmed/28007376"
gene_source <- "Wang_et_al_2017"

# Description of the data:
data_description <- c(
  "977 epilepsy-associated genes were compiled from",
  "the literature and online databases. Genes were",
  "classified into 4 phenotypic categories.",
  "(1) Syndromic epilepsy genes (84 genes).",
  "(2) Neurodevelopment-associated epilepsy",
  "genes (73 genes). (3) Epilepsy-related genes",
  "associated with both physical or other systemic",
  "abnormalities and epilepsy or seizures (536 genes).",
  "(4) Putative epilepsy genes (284 genes)."
)

# Load functions.
devtools::load_all()

# Other imports.
suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
  library(dplyr)
  library(getPPIs)
})

# Directories.
root <- getrd()
gmtdir <- file.path(root, "gmt")
funcdir <- file.path(root, "R")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "downloads")

# Load the data.
myfile <- file.path(downdir, "Wang_et_al_2017_Epilepsy_Genes.xlsx")
data <- read_excel(myfile)

# Clean up the columns.
colnames(data) <- c("Epilepsy", "Neurodevelopment", "Epilepsy-Related", "Potential", "All")

# Melt.
df <- reshape2::melt(data, measure.vars = colnames(data))
colnames(df) <- c("Category", "Gene")

# Map human genes to entrez.
# Why is mapIds returning a list??
genes <- unique(df$Gene)
entrez <- getIDs(genes, from = "symbol", to = "entrez", species = "human")
df$Entrez <- entrez[df$Gene]

# Remove NA.
df <- df[!is.na(df$Entrez), ]

# Map human genes to mouse genes.
df$msEntrez <- getHomologs(df$Entrez, taxid = 10090)

# Remove NA.
df <- df[!is.na(df$msEntrez), ]

# Save table.
myfile <- file.path(tabsdir, paste0(script, ".csv"))
fwrite(df, myfile)

# Split into categories.
gene_groups <- split(df$msEntrez, df$Category)

# Progress report:
ngenes <- length(unique(df$msEntrez))
message(paste("Compiled", ngenes, "mouse epilepsy genes."))
sizes <- sapply(gene_groups, length)
knitr::kable(t(sizes), row.names = FALSE, format = "markdown")

# Save as gmt.
myfile <- file.path(gmtdir, script)
write_gmt(gene_groups, gmt_source = paper_url, gmt_file = myfile)

# Save as rda and generate documentation.
documentDataset(myfile, short_name, Rdir = file.path(root, "R"), datadir)
