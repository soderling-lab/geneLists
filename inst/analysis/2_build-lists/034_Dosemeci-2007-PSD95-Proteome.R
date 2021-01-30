#!/usr/bin/env Rscript

## Parameters:
short_name <- "dosemeci2007PSD95proteome"
script <- "034_Dosemeci-2007-PSD95-Proteome"
ref_url <- "http://www.ncbi.nlm.nih.gov/pubmed?term=17623647"
data_url <- c(
  S1 = "https://www.mcponline.org/highwire/filestream/33409/field_highwire_adjunct_files/0/AffinPur_Top50_SpTr_SuppTable1.xls",
  S2 = "https://www.mcponline.org/highwire/filestream/33409/field_highwire_adjunct_files/4/ParentTop50_SpTr_SuppleTable2.xls",
  S3 = "https://www.mcponline.org/highwire/filestream/33409/field_highwire_adjunct_files/5/Parent_Proteins_SupplementalTable3.xls",
  S4 = "https://www.mcponline.org/highwire/filestream/33409/field_highwire_adjunct_files/3/ParentPeptides_SupplementalTable4.xls",
  S5 = "https://www.mcponline.org/highwire/filestream/33409/field_highwire_adjunct_files/2/AffinPurif_Proteins_SupplementalTable5.xls",
  S6 = "https://www.mcponline.org/highwire/filestream/33409/field_highwire_adjunct_files/1/AffinPurif_Peptides_SupplementalTable6.xls"
)

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

# Download the data.
myfiles <- file.path(downdir, basename(data_url))
download.files(data_url, destdir = downdir, quiet = TRUE)

# Load the data.
all_data <- lapply(myfiles, function(x) read_excel(x, skip = 3))
names(all_data) <- tools::file_path_sans_ext(basename(myfiles))

# Lets focus on T1 and T2.
data <- all_data[c(1, 2)]

# Loop to map ids.
for (i in c(1:2)) {
  ids <- data[[i]]$Entry_ID
  IDs <- file.path(getwd(), "IDs.txt")
  writeLines(ids, IDs)
  # Call python script to map identifiers.
  pyScript <- file.path(root, "Py", "mapIds.py")
  cmd <- paste(pyScript, IDs, "ACC", "P_ENTREZGENEID", "--nsteps 1")
  result <- system(cmd, intern = TRUE)
  # Remove temporary file.
  unlink(IDs)
  # Parse the result.
  entrez <- strsplit(gsub("\\[|\\]|'", "", result), ", ")[[1]]
  entrez[entrez == "None"] <- NA
  # Status.
  not_mapped <- 100 * sum(is.na(entrez)) / length(entrez)
  message(paste("Percent mapped genes:", round(100 - not_mapped, 3), "%."))
  # Map to mouse.
  msEntrez <- getHomologs(entrez, species = "mouse")
  data[[i]] <- tibble::add_column(data[[i]], entrez, .after = "Entry_ID")
  data[[i]] <- tibble::add_column(data[[i]], msEntrez, .after = "entrez")
  # Drop NA.
  data[[i]] <- data[[i]][!is.na(data[[i]]$msEntrez), ]
}

# New names.
names(data) <- sapply(data, function(x) unique(x$Preparation))

# Collect list of genes.
# They are equivalent!
gene_list <- lapply(data, function(x) unique(x$msEntrez))
gene_list[["All"]] <- unique(c(data[[1]]$msEntrez, data[[2]]$msEntrez))

# Just keep one.
gene_list <- list(gene_list[[1]])
names(gene_list) <- "Affinity Purified PSD95-Complex"
message(paste(
  length(gene_list[[1]]), "proteins identified from affinity",
  "purified PSD95 protein complexes."
))

# Save as gene list.
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, ref_url, myfile)

# Save as rda and generate documentation.
documentDataset(myfile, short_name, Rdir = file.path(root, "R"), datadir)
