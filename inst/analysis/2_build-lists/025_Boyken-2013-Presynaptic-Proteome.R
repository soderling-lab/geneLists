#!/usr/bin/env Rscript

## Scrape Synaptic vesicle proteome from Takamori et al., 2006

## Parameters:
short_name <- "boyken2013presynapse"
script <- "025_Boyken-2013-Presynaptic-Proteome"
ref_url <- "https://www.ncbi.nlm.nih.gov/pubmed/23622064"
data_url <- c(
  S1 = "https://ars.els-cdn.com/content/image/1-s2.0-S0896627313001852-mmc2.xlsx",
  S5 = "https://ars.els-cdn.com/content/image/1-s2.0-S0896627313001852-mmc4.xlsx"
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
download.file(data_url[1], myfiles[1])
download.file(data_url[2], myfiles[2])
data <- list()
data[["S1"]] <- read_excel(myfiles[1])
data[["S5"]] <- read_excel(myfiles[2])

## S1 contains data from docked and free synaptic vesicles.
S1 <- data$S1

# Get docked vesicles.
idx <- S1$"ratio docked SV/SV" == "docked SV"
idx[is.na(idx)] <- FALSE
S1 <- S1[idx, ]

# Use gi to map missing ids to entrez.
gi <- sapply(strsplit(S1$"GI number", "gi\\|"), "[", 2)
IDs <- file.path(getwd(), "IDs.txt")
writeLines(gi, IDs)

# Call python script to map identifiers.
pyScript <- file.path(root, "Py", "mapIds.py")
cmd <- paste(pyScript, IDs, "P_GI", "P_ENTREZGENEID")
result <- system(cmd, intern = TRUE)

# Remove temporary file.
unlink(IDs)

# Parse the result.
entrez <- strsplit(gsub("\\[|\\]|'", "", result), ", ")[[1]]
entrez[entrez == "None"] <- NA

# Status.
not_mapped <- 100 * sum(is.na(entrez)) / length(entrez)
message(paste("Percent mapped genes:", round(100 - not_mapped, 3), "%."))

# Map rat entrez to mouse.
msEntrez <- getHomologs(entrez, species = "mouse")

# Add to data.
S1$entrez <- entrez
S1$msEntrez <- msEntrez

# Drop NA.
S1 <- S1[!is.na(S1$msEntrez), ]

## S5 is VGLUT and VGAT data.
S5 <- data$S5
S5 <- filter(S5, S5$"protein name" != "NA")

# Assign vesicles as VGLUT or VGAT.
S5$"average enrichment" <- as.numeric(S5$"average enrichment")
S5$SV_Type <- "VGLUT"
S5$SV_Type[S5$"average enrichment" < 1] <- "VGAT"


# Use gi to map missing ids to entrez.
gi <- sapply(strsplit(S5$"GI number", "gi\\|"), "[", 2)
IDs <- file.path(getwd(), "IDs.txt")
writeLines(gi, IDs)

# Call python script to map identifiers.
pyScript <- file.path(root, "Py", "mapIds.py")
cmd <- paste(pyScript, IDs, "P_GI", "P_ENTREZGENEID")
result <- system(cmd, intern = TRUE)

# Remove temporary file.
unlink(IDs)

# Parse the result.
entrez <- strsplit(gsub("\\[|\\]|'", "", result), ", ")[[1]]
entrez[entrez == "None"] <- NA

# Status.
not_mapped <- 100 * sum(is.na(entrez)) / length(entrez)
message(paste("Percent mapped genes:", round(100 - not_mapped, 3), "%."))

# Map rat entrez to mouse.
msEntrez <- getHomologs(entrez, species = "mouse")

# Add to data.
S5$entrez <- entrez
S5$msEntrez <- msEntrez

# Drop NA.
S5 <- S5[!is.na(S5$msEntrez), ]

# Collect list of genes.
gene_list <- list()
gene_list[["Docked SV"]] <- unique(S1$msEntrez)
gene_list[["VGLUT SV"]] <- unique(S5$msEntrez[S5$SV_Type == "VGLUT"])
gene_list[["VGAT SV"]] <- unique(S5$msEntrez[S5$SV_Type == "VGAT"])

# Summary
df <- data.frame("N Genes" = sapply(gene_list, length))
df <- tibble::add_column(df, "Class" = rownames(df), .before = 1)
knitr::kable(df, row.names = FALSE)

# Save as gene list.
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, ref_url, myfile)

# Save as rda and generate documentation.
documentDataset(myfile, short_name, Rdir = file.path(root, "R"), datadir)
