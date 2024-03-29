#!/usr/bin/env Rscript

## Building a DBD-associated gene collection.

## DisGeneNet Datasets:
# [1] Curated_Disease_Genes
# [2] All_Disease_Genes
# [3] Curated_Variants
# [4] All_Variants

root <- getrd()
renv::load(root, quiet = TRUE)

## Parameters:
script <- "047_Combined-DisGeneNET-NDD-geneSet"
short_name <- "disgeneNDD"
datasets <- c("Curated_Disease_Genes", "All_Disease_Genes", "Curated_Variants", "All_Variants")
diseases <- c(
  "alzheimer",
  "multiple sclerosis",
  "amyotrophic lateral sclerosis",
  "nerve degeneration",
  "lewy body disease",
  "friedreich ataxia",
  "parkinson disease",
  "huntington disease",
  "spinal muscular atrophy"
)

# Imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(getPPIs)
  library(parallel)
  library(doParallel)
  library(data.table)
})

# Load functions in projects R/ directory.
suppressWarnings({
  devtools::load_all()
})

# Register some workers for parallel execution.
nThreads <- detectCores() - 1
workers <- makeCluster(nThreads)
registerDoParallel(workers)

# Directories.
root <- getrd()
gmtdir <- file.path(root, "gmt")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

# Define a function that gets DisGeneNet data:
getDisGeneNETdata <- function(dataset = c(
                                "Curated_Disease_Genes", "All_Disease_Genes",
                                "Curated_Variants", "All_Variants"
                              ),
                              sources = c(
                                "CTD_human", "BEFREE", "PSYGENET", "LHGDN", "HPO",
                                "GENOMICS_ENGLAND", "RGD", "GWASCAT", "GWASDB", "MGD",
                                "CLINVAR", "UNIPROT", "CTD_rat", "ORPHANET", "CTD_mouse"
                              ),
                              base_url = "https://www.disgenet.org/static/disgenet_ap1/files/downloads") {
  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(getPPIs)
  })
  # Download, unzip, and load the DisGeneNet data.
  # Variant Gene map is used for mapping gene variants.
  datasets <- c(
    Curated_Disease_Genes = "curated_gene_disease_associations.tsv.gz",
    All_Disease_Genes = "all_gene_disease_associations.tsv.gz",
    Curated_Variants = "curated_variant_disease_associations.tsv.gz",
    All_Variants = "all_variant_disease_associations.tsv.gz",
    Variant_Gene_Map = "variant_to_gene_mappings.tsv.gz"
  )
  myurl <- file.path(base_url, datasets[dataset])
  download.file(myurl, destfile = basename(myurl), quiet = TRUE)
  system(command = paste("gunzip", basename(myurl)))
  myfile <- tools::file_path_sans_ext(basename(myurl))
  data <- fread(myfile)
  unlink(myfile)
  # If working with gene variants, get mapping data from DisGeneNet.
  if (grepl("Variants", dataset)) {
    myurl <- file.path(base_url, datasets["Variant_Gene_Map"])
    download.file(myurl, destfile = basename(myurl), quiet = TRUE)
    system(command = paste("gunzip", basename(myurl)))
    myfile <- tools::file_path_sans_ext(basename(myurl))
    varmap <- fread(myfile)
    unlink(myfile)
    # Add geneID and geneSymbol columns to data.
    entrez <- varmap$geneId[match(data$snpId, varmap$snpId)]
    genes <- varmap$geneSymbol[match(data$snpId, varmap$snpId)]
    data <- tibble::add_column(data, geneId = entrez, .after = 1)
    data <- tibble::add_column(data, geneSymbol = genes, .after = 2)
    data <- data[geneId != "NA"] # Remove unmapped variants.
  }
  # DisGeneNET is a database of human genes and their associated diseases.
  # Map human genes in to their mouse homologs.
  hsEntrez <- data$geneId
  msEntrez <- getHomologs(hsEntrez, species = "mouse")
  nMsGenes <- length(unique(msEntrez))
  data <- tibble::add_column(data, msEntrez = msEntrez, .after = 1)
  # Remove rows with unmapped genes.
  data <- subset(data, !is.na(data$msEntrez))
  # Split rows (genes) derived from multiple sources.
  data <- data %>% tidyr::separate_rows(source, sep = ";")
  # Keep the data of interest.
  data <- subset(data, data$source %in% sources)
  return(data)
}

# Download and clean-up the data.
data <- lapply(datasets, function(x) getDisGeneNETdata(x))
names(data) <- datasets

# Combine datasets by shared column names.
colNames <- Reduce(intersect, sapply(data, colnames))

# Note: The following columns will be dropped:
# allNames <- Reduce(union,sapply(data,colnames))
# allNames[allNames %notin% colNames]

# Merge dataframes.
# This takes a couple seconds...
data <- data %>% purrr::reduce(full_join, by = colNames)

# If data was found in multiple dataframes then a column is added.
# e.g. chromosome.x; remove these duplicated columns.
data <- data[, -c(grep("\\.", colnames(data)))]

# Save the data.
myfile <- file.path(tabsdir, "DisGENEnet.csv")
fwrite(data, myfile)

# Write a function to find rows for diseases of interest.
getRows <- function(diseases, i) {
  idx <- grep(diseases[i], tolower(data$diseaseName))
  return(idx)
}

# Parallelize the task with dopar.
rows <- foreach(i = seq_along(diseases)) %dopar% {
  getRows(diseases, i)
}

# Send workers home.
suppressWarnings(stopCluster(workers))

# Split data into disease groups.
data_list <- lapply(rows, function(idx) data[idx, ])
names(data_list) <- diseases

# Check disorder group sizes.
sizes <- sapply(data_list, function(x) length(unique(x$msEntrez)))
message("\nSummary of gene-disease associations:")
knitr::kable(t(sizes), row.names = FALSE, format = "markdown")
message("\n")

# Write as gmt file.
gmt_list <- lapply(data_list, function(x) x$msEntrez)
gmt_file <- file.path(gmtdir, paste0(script, ".gmt"))
write_gmt(gmt_list, gmt_source = "DisGeneNET", gmt_file)

# Status report.
nGenes <- length(unique(unlist(sapply(data_list, function(x) x$msEntrez))))
nDisorders <- length(diseases)
message(paste0(
  "Compiled ", nGenes, " mouse genes associated with ",
  nDisorders, " Neurodegenerative disorders!"
))

# Save as rda and generate documentation.
documentDataset(gmt_file, short_name, Rdir = file.path(root, "R"), datadir)
