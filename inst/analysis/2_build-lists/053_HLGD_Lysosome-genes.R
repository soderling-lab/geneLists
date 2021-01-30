#!/usr/bin/env Rscript
# Compile lysosome genes/proteins from the (Human and Mouse) Lysosome Gene
# Database: http://lysosome.unipg.it/

## Parameters:
script <- "053_HGLD_Lysosome-genes"
short_name <- "hlgd"
# NOTE: data were manually downloaded as csv. 1 Error in mus data was fixed manually.
data_files <- c("mouse_LGD.csv", "human_LGD.csv") # in root/downloads
ignore <- c(10802651, 15044231, 10592173, 17367534) # ignore Uniprot, GO, KEGG, Reactome  datasets--indirect data sources.

# Load renv.
root <- dirname(dirname(getwd()))

renv::load(root, quiet = TRUE)

# Imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(getPPIs)
  library(data.table)
})

# Load functions in projects R/ directory.
suppressWarnings({ devtools::load_all() })

# Project directories.
gmtdir <- file.path(root, "gmt")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "downloads")

# Define a function that loads the manually downloaded data in root/downloads.
load_data <- function(data_file) {
  data_file <- file.path(downdir, data_files[2])
  # Load the data from LGD. 
  # To make sure the formatting is consistent I did it this round-about way
  # using readLines.
  raw_data <- readLines(data_file)[-1]
  raw_list <- strsplit(raw_data, "\t")
  idx <- which(sapply(raw_list, length) != 4)
  to_fix <- split(idx, rep(seq(1, length(idx) / 2), each = 2))
  fixed <- lapply(to_fix, function(x) {
    idx <- x[1]
    idy <- x[2]
    fixed <- c(raw_list[[idx]], raw_list[[idy]])
    return(fixed[fixed != ""])
  })
  fixed_list <- c(raw_list[-idx], fixed)
  df <- as.data.frame(do.call(rbind, fixed_list))
  colnames(df) <- c("Entrez", "Gene", "Name", "Ref")
  rownames(df) <- NULL
  dt <- as.data.table(tidyr::separate_rows(df, "Ref", sep = ","))
  return(dt)
}

# Warnings about incomplete final lines can be ignored.
suppressWarnings({
	data_list <- lapply(file.path(downdir, data_files), load_data)
})
names(data_list) <- sapply(strsplit(data_files, "_"), "[", 1)

# Combine the two datasets.
data_dt <- bind_rows(data_list, .id = "Species")
subdt <- data_dt %>% filter(Ref %notin% ignore)

# Check the refs.
refs_dt <- data.table(ref = unique(subdt$Ref))
ref_list <- list( # Manually inspect each reference:
  "20957757" = "The proteome of lysosomes",
  "19556463" = "A gene network regulating lysosomal biogenesis and function",
  "18977398" = "Proteomics of the lysosome",
  "17897319" = "Integral and associated lysosomal membrane proteins",
  "16709564" = "Identification and validation of mannose 6-phosphate...",
  "17258946" = "The human urine mannose 6-phosphate glycoproteome",
  "15789345" = "The human brain mannose 6-phosphate glycoproteome...",
  "16399764" = "Identification of sites of mannose 6-phosphorylation...",
  "11079561" = "Towards a human repertoire of monocytic lysosomal proteins",
  "12203898" = "Proteomic analysis of human lysosomes...",
  "9859869" = "Two-dimensional mapping and microsequencing of lysosomal...",
  "21752829" = "Characterization of the CLEAR network reveals..."
)

# Collapse references.
data <- subdt %>%
  group_by(Species, Entrez) %>%
  summarize(Refs = paste(Ref, collapse = ";"))

# Extract mouse data.
ms_dt <- data %>% filter(Species == "mouse")

# Extract human data, and then
# map human entrez to mouse homologs, this takes several moments.
hs_dt <- data %>% filter(Species == "human")
msEntrez <- getPPIs::getHomologs(hs_dt$Entrez, species = "mouse")
hs_dt$msEntrez <- msEntrez

# Drop unmapped genes from human subdat.
hs_dt <- hs_dt %>%
  filter(!is.na(msEntrez)) %>%
  select(Species, msEntrez, Refs)

# Fix column names.
colnames(hs_dt)[2] <- "Entrez"
# Insure that Entrez is character.
hs_dt$Entrez <- as.character(hs_dt$Entrez)

# Combine with mouse data.
gene_dt <- bind_rows(list(ms_dt, hs_dt)) %>% filter(!is.na(Entrez))

# Coerce to gene list format.
gene_list <- list("LGD" = gene_dt$Entrez)

# Status report.
nGenes <- length(unique(gene_dt$Entrez))
nPub <- length(refs_dt$ref)
message(paste("Compiled", nGenes, "mouse Lysosome genes from", nPub, "studies."))

# Write as gmt file.
gmt_file <- file.path(gmtdir, paste0(script, ".gmt"))
write_gmt(gene_list, gmt_source = "LGD", gmt_file)

# Save as rda and generate documentation.
documentDataset(gmt_file, short_name, Rdir = file.path(root, "R"), datadir)
