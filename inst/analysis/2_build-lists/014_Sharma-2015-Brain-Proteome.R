#!/usr/bin/env Rscript

## Scrape mouse brain proteome data from Sharma et al., 2015.
# Note: I can't seem replicate their reported number of genes that are
# enriched in the brain or brain structures.

renv::load(getrd())

## Parameters:
script <- "014_Sharma-2015-Brain-Proteome"
short_name <- "sharma2015brain"
url <- "http://141.61.102.17/mousebrainproteome/datasets/BrainRegions/BrainRegionTotalDataset_Log2FoldChange.xls"
paper <- "https://www.ncbi.nlm.nih.gov/pubmed/26523646"

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
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "downloads")

# Download the data.
myfile <- file.path(downdir, basename(url))
download.file(url, myfile, quiet = TRUE)

# Change file extension.
# Although readxl should be able to handle .xls, trying to read the
# file as is with read_excel or read_xls results in the error:
# libxls error: Unable to open file libxls -- probably need the libxls library.
# Changing the extension seems to work fine.
newfile <- paste0(tools::file_path_sans_ext(myfile), ".xlsx")
invisible(file.copy(myfile, newfile))

# Load the data.
# Changing file extension might result in loss of second sheet, but
# we can recover the DA proteins.
data <- read_excel_sheets(newfile)
data <- data[[1]]

# Clean up the data.
col_names <- paste0(as.character(data[1, ]), seq(ncol(data)))
colnames(data) <- col_names
data <- data[2:nrow(data), ]

# Split rows with multiple gene names.
data <- data %>% tidyr::separate_rows("Gene names1", sep = ";")

# Split rows with multiple protein ids.
data <- data %>% tidyr::separate_rows("Majority protein IDs58", sep = ";")

# Map uniprot to entrez.
uniprot <- data$"Majority protein IDs58"
uniprot <- sapply(strsplit(uniprot, "-"), "[", 1) # Ignore isoforms.
entrez <- getIDs(uniprot, from = "uniprot", to = "entrez", species = "mouse")

# Map missing entrez using gene symbols.
# Why does mapIds return a list?
is_missing <- is.na(entrez)
genes <- data$"Gene names1"[is_missing]
entrez[is_missing] <- getIDs(genes, from = "symbol", to = "entrez", species = "mouse")

# What percent of ids are missing?
is_missing <- is.na(entrez)
# 100 * (sum(is_missing) / length(entrez)) # 1.5% Not terrible...

# Try mapping remaining missing uniprot ids.
myfile <- file.path(downdir, "not_mapped.csv")
writeLines(uniprot[is_missing], myfile)

# Map ids with mgi batch tool.
myfile <- file.path(downdir, "MGIBatchReport_20200301_113813.xlsx")
idmap <- read_excel(myfile)

# Collect mgi ids.
idy <- grep("MGI", colnames(idmap))
idx <- match(uniprot[is_missing], idmap$Input)
mgi <- idmap[[idy]][idx]
names(mgi) <- uniprot[is_missing]
mgi[grep("No associated gene", mgi)] <- NA
entrez[is_missing] <- getIDs(mgi, from = "mgi", to = "entrez", species = "mouse")

# Add entrez to data.
data <- tibble::add_column(data, entrez, .after = "Gene names1")

# Remove unmapped genes.
data <- data %>% filter(!is.na(entrez))

# All genes.
all_genes <- unique(data$entrez)
N <- length(all_genes)

# Get significantly differentially abundant proteins.
cols <- c(
  "Brainstem24", "Cerebellum25", "Corpus Callosum26", "Motor Cortex27",
  "Olfactory Bulb28", "Optic Nerve29", "Prefrontal Cortex30", "Striatum31",
  "Thalamus32", "Hippocampus33"
)
idy <- c(
  grep("entrez", colnames(data)),
  grep(">4 fold", colnames(data)),
  match(cols, colnames(data))
)
subdat <- data[, idy]

# Clean-up a little.
colnames(subdat)[2] <- "Sig_Brain"

# Fix significance column.
subdat$Sig_Brain <- !is.na(subdat$Sig_Brain)

# Get significant proteins.
# Melt into tidy df.
df <- subdat %>%
  filter(Sig_Brain) %>%
  reshape2::melt(id.var = "entrez")

# Fix significance column.
df$value <- !is.na(df$value)

# Fix-up column names.
colnames(df) <- c("entrez", "region", "sig")

# Remove numbers from brain region names.
df$region <- gsub("[0-9]", "", df$region)

# Sig_brain indicates that the protein is enriched in any brain region.
df$region[df$region == "Sig_Brain"] <- "Brain"

# Drop any ns proteins.
df <- df %>% filter(df$sig)

# Split into lists.
data_list <- df %>%
  group_by(region) %>%
  group_split()
names(data_list) <- sapply(data_list, function(x) unique(x$region))

# Remove duplicates.
gene_list <- lapply(data_list, function(x) unique(x$entrez))

# Summary.
df <- data.frame("Genes" = sapply(gene_list, length))
df <- tibble::add_column(df, "Brain Region" = rownames(df), .before = 1)
knitr::kable(df, row.names = FALSE)

# Write as gmt.
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, gmt_source = paper, myfile)

# Save as rda and generate documentation.
documentDataset(myfile, short_name, Rdir = file.path(root, "R"), datadir)
