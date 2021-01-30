#!/usr/bin/env Rscript

## Scrape data from Velmeshev et al., 2019.

renv::load(getrd())

## Parameters
script <- "011_Velmeshev-2019-ASD-geneSet"
short_name <- "velmeshev2019ASD"
paper <- "https://www.ncbi.nlm.nih.gov/pubmed/31097668"

# Imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(getPPIs)
  library(data.table)
  library(org.Hs.eg.db)
})

# Functions.
suppressWarnings({
  invisible({
    devtools::load_all()
  })
})

# Directories.
root <- getrd()
gmtdir <- file.path(root, "gmt")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "downloads")

# Urls to supplemental data.
urls <- c(
  S1 = "https://science.sciencemag.org/highwire/filestream/726807/field_highwire_adjunct_files/0/aav8130_Data-S1.xlsx",
  S2 = "https://science.sciencemag.org/highwire/filestream/726807/field_highwire_adjunct_files/1/aav8130_Data-S2.xlsx",
  S3 = "https://science.sciencemag.org/highwire/filestream/726807/field_highwire_adjunct_files/2/aav8130_Data-S3.xls",
  S4 = "https://science.sciencemag.org/highwire/filestream/726807/field_highwire_adjunct_files/3/aav8130_Data-S4.xls",
  S5 = "https://science.sciencemag.org/highwire/filestream/726807/field_highwire_adjunct_files/4/aav8130_Data-S5.xlsx"
)

# Description of the data.
descriptions <- list(
  S1 = "Sample and clinical information for ASD and epilepsy individuals.",
  S2 = "List of captured nuclei and associated metadata for ASD and epilepsy cohorts.",
  S3 = "List of cluster-specific and regional gene markers.",
  S4 = c(
    "List of cell type-specific genes differentially expressed in ASD and epilepsy,",
    " as well as region-specific and individual-specific gene expression changes in ASD."
  ),
  S5 = c(
    "Results of whole-exome sequencing analysis of ASD patients.",
    " The first tab includes high-confidence variants, the second tab",
    " includes variants that are associated with downregulation of",
    " corresponding genes in the same ASD patient when compared to control samples;",
    " the other tabs include unfiltered lists of variants for each individual with less",
    " confidence in association with ASD, epilepsy or psychiatric disease."
  )
)

# Download the data.
myfiles <- basename(urls)
download.files(urls, quiet = TRUE)

#---------------------------------------------------------------------
## Load the data.
#---------------------------------------------------------------------

## Sheet 1 of S1 has unique formatting. Load it first.
myfile <- list.files(pattern = "Data-S1")
sheets <- excel_sheets(myfile)
S1 <- list()
colNames <- readLines(file.path(downdir, "colnames.txt"))
S1[[1]] <- read_excel(myfile, sheet = 1, skip = 2, col_names = colNames)
S1 <- c(S1, lapply(sheets[-1], function(x) read_excel(myfile, sheet = x)))
names(S1) <- sheets

# Load the rest of the data. Combine with S1 into a list.
data <- list()
data[[1]] <- S1
data <- c(data, lapply(myfiles[-1], read_excel_sheets))
names(data) <- names(urls)

# Remove downloaded files.
unlink(basename(urls))

# Table S4 contains DE genes that are:
# Cell-type specific.
# Region specific (ACC and PFC).
# WGCNA modules.
# Cell-type specific from specific individuals.

# Lets just collect all the genes and see what they are.
S4 <- data$S4
df <- S4[["ASD_DEGs"]]

# Map ensembl ids.
db <- org.Hs.egENSEMBL
ensembl_ids <- unlist(as.list(db[mappedkeys(db)]))
idx <- match(df$"gene ID", ensembl_ids)
entrez <- names(ensembl_ids[idx])
df <- tibble::add_column(df, entrez, .after = "gene ID")

# Remove unmapped ids.
out <- is.na(df$entrez)
df <- df[!out, ]
message(paste("Number of genes that were not mapped:", sum(out)))

# Map to mouse genes.
msEntrez <- getHomologs(df$entrez, taxid = 10090)
df <- tibble::add_column(df, msEntrez, .after = "entrez")

# Remove unmapped ids.
out <- is.na(df$msEntrez)
df <- df[!out, ]

# Save table.
myfile <- file.path(tabsdir, paste0(script, ".csv"))
fwrite(df, myfile)

# Summary.
nGenes <- length(unique(df$msEntrez))
message(paste(
  "Compiled", nGenes, "mouse genes that are differentially",
  "expressed in humans with ASD."
))

# Collect genes in a list.
gene_list <- split(df$msEntrez, df$"Cell type")
names(gene_list) <- paste(names(gene_list), "DEGs")
gene_list[["ASD DEGs"]] <- unique(df$msEntrez)
idx <- df$"Fold change" > 0
gene_list[["ASD Up-regulated DEGs"]] <- df$msEntrez[idx]
gene_list[["ASD Down-regulated DEGs"]] <- df$msEntrez[!idx]

# Status report.
sizes <- sapply(gene_list, length)
mytable <- data.table("N Genes" = sizes)
mytable <- tibble::add_column(mytable, "Class" = names(sizes), .before = 1)
knitr::kable(mytable, row.names = FALSE)

# Save as gmt.
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, paper, myfile)

# Save as rda and generate documentation.
documentDataset(myfile, short_name, Rdir = file.path(root, "R"), datadir)
