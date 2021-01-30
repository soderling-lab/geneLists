#!/usr/bin/env Rscript

#' ---
#' title: Compilation of CORUM protein complexes.
#' description:
#' authors: Tyler W Bradshaw
#' ---

## Compile CORUM protein complexes.
identity_threshold <- 0.9

## Parameters:
script <- "043_CORUM-Complexes"
short_name <- "corum"
tags <- c("corum", "protein complex", "complexes", "database")
ref_url <- "https://www.ncbi.nlm.nih.gov/pubmed/30357367"
data_url <- "https://mips.helmholtz-muenchen.de/corum/download/coreComplexes.txt.zip"

# Load renv.
renv::load(getrd())

# Imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(getPPIs)
  library(TBmiscr)
  library(data.table)
})

# Load functions in root/R.
load_all()

# Directories.
root <- getrd()
gmtdir <- file.path(root, "gmt")
datadir <- file.path(root, "data")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "downloads")

# Download Chorum complexes.
zip_file <- file.path(downdir, basename(data_url))
myfile <- tools::file_path_sans_ext(zip_file)
message(paste("Downloading PPI Complexes from CORUM database..."))
download.file(data_url, zip_file, quiet = TRUE)

# Extract zipped file and remove temporary files.
unzip(zip_file, exdir = downdir)
data <- data.table::fread(myfile)
unlink(zip_file)
unlink(myfile)

## Clean-up the Corum data.
# We need to extract Entrez IDs associated with protein complexes.
# Split Entrez ID column.
data <- tidyr::separate_rows(data, "subunits(Entrez IDs)", sep = ";")

# Let's collect the relevant data columns.
cols <- c(
  "ComplexID", "subunits(Entrez IDs)", "Organism", "PubMed ID",
  "ComplexName", "Synonyms",
  "Protein complex purification method", "Complex comment"
)
corum <- data %>% select(all_of(cols))

# Clean-up colnames.
colnames(corum) <- c(
  "ComplexID", "Entrez", "Organism", "PMID", "Name",
  "Synonym", "Method", "Comment"
)

# Map entrez ids to mouse homologs.
musEntrez <- getHomologs(corum$Entrez, species = "mouse")

# Add mouse homologs to data.
corum <- tibble::add_column(corum, musEntrez, .after = "Entrez")

# Summarize total number of proteins per complex as well as the number of mouse
# homologs per complex and PMIDs.
corum_complexes <- corum %>%
  group_by(ComplexID) %>%
  summarize(
    Name = unique(Name),
    nProts = length(unique(Entrez)),
    nHomologs = sum(!is.na(unique(musEntrez))),
    PMIDs = paste(unique(PMID), collapse = "; "),
    Entrez = paste(unique(Entrez), collapse = "; "),
    musEntrez = paste(unique(musEntrez[!is.na(musEntrez)]),
      collapse = "; "
    )
  )

# Separate out mouse Entrez column again.
corum_complexes <- tidyr::separate_rows(corum_complexes, "musEntrez", sep = "; ")

# Add percent coverage.
p <- corum_complexes$nHomologs / corum_complexes$nProts
corum_complexes <- tibble::add_column(corum_complexes,
  percentID = p, .after = "nHomologs"
)

# Subset the data: complete identity.
subdat <- corum_complexes %>% filter(percentID >= identity_threshold)
n <- length(unique(subdat$Name))
message(paste("Collected", formatC(n, big.mark = ","), "mouse complexes."))

# Split into gene list.
data_list <- subdat %>%
  group_by(Name) %>%
  group_split()

# Collect list of genes.
gene_list <- sapply(data_list, function(x) unique(x$musEntrez))
names(gene_list) <- sapply(data_list, function(x) unique(x$Name))

# Save as gmt, and then save as rda and generate documentation.
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, ref_url, myfile)
documentDataset(myfile, short_name, Rdir = file.path(root, "R"), datadir)
