#!/usr/bin/env Rscript

## Scrape Synaptic proteins from the SynSysNet database.
# Website is scraped with python script.

## Parameters:
short_name <- "synsysnet"
script <- "021_SynSysNet-Proteome"
ref_url <- "http://bioinformatics.charite.de/synsys/index.php?site=syn_class"

# Load renv:
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
gmtdir <- file.path(root, "gmt")
datadir <- file.path(root, "data")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "downloads")

# Load the data.
myfile <- file.path(downdir, "SynSysNet-Synaptic-Protein-Tree.csv")
data <- fread(myfile)

# Collect genes.
data_list <- split(data$V2, data$V1)
data_list <- lapply(data_list, function(x) unlist(strsplit(x, ",")))

# Drop the first element, its empty.
data_list[[1]] <- NULL

# Collect all pre,post, and unspecified proteins.
data_list[["post"]] <- unlist(data_list[grep("post", names(data_list))], use.names = FALSE)
data_list[["pre"]] <- unlist(data_list[grep("pre", names(data_list))], use.names = FALSE)
data_list[["unspec"]] <- unlist(data_list[grep("unspec", names(data_list))], use.names = FALSE)

# Map Uniprot IDs to entrez.
entrez <- lapply(data_list, function(x) {
  getIDs(x, from = "uniprot", to = "entrez", species = "human")
})

# Attrition is pretty bad. Try mapping directly with uniprot db.
uniprot <- unique(unlist(data_list))
writeLines(uniprot, file.path(downdir, paste0(Sys.Date(), "_uniprot.csv")))

# Load mapping table from uniprot.
myfile <- file.path(downdir, "2020-03-03_uniprot-map.xlsx")
protmap <- read_excel(myfile)
colnames(protmap)[1] <- "Uniprot"

# Remove rows with missing values as these will rais errors.
out <- is.na(protmap$Status)
protmap <- protmap[!out, ]

# Separate rows with multiple Uniprot ids.
protmap <- protmap %>% tidyr::separate_rows(Uniprot, sep = ",")

# Separate rows with multiple gene names.
protmap <- protmap %>% tidyr::separate_rows("Gene names", sep = " ")

# Map ids.
entrez <- getIDs(protmap$Entry, from = "uniprot", to = "entrez", species = "human")

# Map missing with gene names.
is_missing <- is.na(entrez)
entrez[is_missing] <- getIDs(protmap$"Gene names"[is_missing],
  from = "symbol", to = "entrez", species = "human"
)

# Try one last time.
is_missing <- is.na(entrez)
entrez[is_missing] <- getIDs(protmap$"Uniprot"[is_missing],
  from = "uniprot", to = "entrez", species = "human"
)

# Check missing.
# Better.
is_missing <- is.na(entrez)
percent_missing <- 100 * (sum(is_missing) / length(entrez))
message(paste("Percent unmapped IDs:", round(percent_missing, 3)))

# Vector of mapped genes.
entrez <- entrez[!is_missing]
names(entrez) <- protmap$Uniprot[!is_missing]

# Loop to map uniprot ids to entrez.
gene_list <- lapply(data_list, function(x) unique(entrez[x]))

# Remove na.
gene_list <- lapply(gene_list, function(x) x[!is.na(x)])

# Map human entrez to mouse.
## FIXME: getHomologs quiet doese not supress output!
gene_list <- lapply(gene_list, function(x) getHomologs(x, species = "mouse", quiet = TRUE))

# Remove na.
gene_list <- lapply(gene_list, function(x) x[!is.na(x)])

# Summary.
df <- data.frame("N Genes" = sapply(gene_list, length))
df <- tibble::add_column(df, "Category" = rownames(df), .before = 1)
knitr::kable(df, row.names = FALSE)

# Save as gene list.
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, ref_url, myfile)

# Save as rda and generate documentation.
documentDataset(myfile, short_name, Rdir = file.path(root, "R"), datadir)
