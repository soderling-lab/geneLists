#!/usr/bin/env Rscript

## Scrape cell type marker gene data from Velmeshev et al., 2019

renv::load(getrd())

## Parameters:
script <- "012_Velmeshev-2019-CellTypes"
short_name <- "velmeshev2019cells"
paper <- "https://www.ncbi.nlm.nih.gov/pubmed/31097668"
url <- "https://science.sciencemag.org/highwire/filestream/726807/field_highwire_adjunct_files/2/aav8130_Data-S3.xls"

# Imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(getPPIs)
  library(data.table)
})

# Examine the number of genes assigned to each cell type.
devtools::load_all()

# Directories.
root <- getrd()
gmtdir <- file.path(root, "gmt")
datadir <- file.path(root, "data")

# Download the data.
download.file(url, destfile = basename(url), quiet = TRUE)

# Load into R.
data <- read_excel(basename(url))

# Remove temporary file.
unlink(basename(url))

# Map human genes to entrez ids.
# FIXME: if ensembl is provided it matches multiple identifier types!
entrez <- getIDs(data$"Gene ID", from = "ensembl$", to = "entrez", species = "human")
not_mapped <- sum(is.na(entrez))
warning(paste(not_mapped, "human ensembl IDs were not successfully mapped."))

# Add to data and drop unmapped genes.
data$Entrez <- entrez
subdat <- data %>% filter(!is.na(Entrez))

# map human entrez to mouse genes.
msEntrez <- getHomologs(subdat$Entrez, species = "mouse")
not_mapped <- sum(is.na(msEntrez))
warning(paste(not_mapped, "human genes do not have a homologous mouse gene."))

# Add to data and drop unmapped genes.
subdat$msEntrez <- msEntrez
filt_data <- subdat %>% filter(!is.na(msEntrez))

# Split into cell type groups.
colnames(filt_data)[1] <- "Cell_type"
data_list <- filt_data %>%
  group_by(Cell_type) %>%
  group_split()
names(data_list) <- sapply(data_list, function(x) unique(x$Cell_type))

# Collect the genes.
gene_list <- sapply(data_list, function(x) x$msEntrez)
names(gene_list) <- names(data_list)

# Summary.
df <- data.frame("Marker_Genes" = sapply(gene_list, length))
df <- tibble::add_column(df, "Cell_Type" = rownames(df), .before = 1)
knitr::kable(df, row.names = FALSE)

# Save gene list as gmt.
gmt_file <- file.path(gmtdir, script)
write_gmt(gene_list, paper, gmt_file)

# Generate documentation.
documentDataset(gmt_file, short_name, Rdir = file.path(root, "R"), datadir)
