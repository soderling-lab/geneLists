#!/usr/bin/env Rscript

# Scrape Brain Proteome Modules from Courtland et al., 2020

# Load renv.
here <- getwd()
root <- dirname(dirname(here))
renv::load(root)

# Urls to data.
data_source <- "https://www.biorxiv.org/content/10.1101/2020.08.06.239517v1"
data_url <- "https://github.com/twesleyb/SwipProteomics/raw/master/data/partition.rda"
map_url <- "https://github.com/twesleyb/SwipProteomics/raw/master/data/gene_map.rda"

script <- "060_Courtland-2020-Brain-Spatial-Proteome"
short_name <- "courtland2020"

# Imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(getPPIs)
  library(data.table)
})

# Functions.
suppressWarnings({
  devtools::load_all()
})

# project directories
gmtdir <- file.path(root, "gmt")
datadir <- file.path(root, "data")

# download the data
myfile <- basename(data_url)
download.file(data_url, myfile, quiet = TRUE)

# load the data
load(myfile) # partition

# remove temporary file
invisible(unlink(myfile))

# download the gene map
myfile <- basename(map_url)
download.file(map_url, myfile, quiet = TRUE)

# load the data
load(myfile) # gene_map

# remove temporary file
invisible(unlink(myfile))

# use gene_map to map uniprot to entrez
idx <- match(names(partition),gene_map$uniprot)
entrez <- gene_map$entrez[idx]

if (sum(is.na(entrez)) != 0) { stop("There are unmapped genes.") }

# split into groups of genes (modules);
# remove M0 - its not a module
modules <- split(entrez,partition)
names(modules) <- paste0("M",names(modules))
modules <- modules[names(modules) != "M0"]

# save as gmt
gene_list <- modules
myfile <- file.path(gmtdir, script)
write_gmt(gene_list, data_source, myfile)

# save as rda and generate documentation
documentDataset(myfile, short_name, Rdir = file.path(root, "R"), datadir)
