#!/usr/bin/env Rscript

## Building a Synaptic GO collection from Koopmans et al., 2019

renv::load(getrd())

## User parameters:
script <- "013_Koopmans-2019-SynGO-geneSet"
short_name <- "synGO"
release <- "20180731" # Most recent SynGO release.
org <- "mouse" # Which organism, human or mouse?
data_source <- "https://www.syngoportal.org/"

# Imports.
suppressPackageStartupMessages({
  library(readxl)
  library(tidyr)
  library(data.table)
  library(dplyr)
  library(getPPIs)
})

# Functions.
suppressWarnings({
  devtools::load_all()
})

# Directories.
root <- getrd()
gmtdir <- file.path(root, "gmt")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "downloads")

# Downoad the raw data.
message("Downloading Synaptic GO (SyngGO) database... \n")
base_url <- "https://syngoportal.org/data/download.php?file=SynGO_bulk_download_release_"
myurl <- paste0(base_url, release, ".zip")
myfile <- file.path(downdir, sapply(strsplit(basename(myurl), "file="), "[", 2))
download.file(myurl, destfile = myfile, method = "curl", quiet = TRUE, extra = "--insecure")
cmd <- paste0("unzip -oq ", myfile, " -d /", tools::file_path_sans_ext(myfile))
system(command = cmd)
unlink(myfile)

# Load SynGO data.
mydir <- tools::file_path_sans_ext(myfile)
myfiles <- list.files(mydir, pattern = "xlsx", full.names = TRUE)
syngo_data <- sapply(myfiles, read_excel)
names(syngo_data) <- gsub(".xlsx", "", basename(names(syngo_data)))

# Get mouse Entrez IDs by mapping MGI IDs to Entrez.
genes_df <- syngo_data$syngo_genes

# Seperate rows with multiple MGI ids.
# Warning caused by missing (NA) MGI ids can be ignored.
genes_df <- genes_df %>% separate_rows(mgi_id, sep = ",")
mgi <- paste0("MGI:", genes_df$mgi_id)
msEntrez <- getIDs(mgi, from = "mgi", to = "entrez", species = "mouse")
genes_df$msEntrez <- msEntrez

# Collect GO annotations and add this column to syngo genes.
anno_df <- syngo_data$syngo_annotations
idx <- match(anno_df$"human ortholog gene hgnc_id", genes_df$hgnc_id)
go_df <- data.frame(
  "id" = as.character(anno_df$"GO term ID"),
  "domain" = as.character(anno_df$"GO domain"),
  "pmid" = as.character(anno_df$PMID),
  "name" = as.character(anno_df$"GO term name"),
  "msEntrez" = as.character(genes_df$msEntrez[idx])
)

# Save to file.
myfile <- file.path(tabsdir, paste0(script, ".csv"))
data.table::fwrite(go_df, myfile)

# Group into gene groups.
data_list <- go_df %>%
  group_by(id) %>%
  group_split()
names(data_list) <- as.character(sapply(data_list, function(x) unique(x$id)))

# Gene list
gene_list <- lapply(data_list, function(x) x$msEntrez)

# Summarize SYNGO terms.
df <- data.frame("Genes" = sapply(gene_list[grep("SYNGO", names(gene_list))], length))
df <- tibble::add_column(df, "SYNGO" = rownames(df), .before = 1)
knitr::kable(df, row.names = FALSE)

# Write gmt.
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, data_source, myfile)

# Save as rda and generate documentation.
documentDataset(myfile, short_name, Rdir = file.path(root, "R"), datadir)
