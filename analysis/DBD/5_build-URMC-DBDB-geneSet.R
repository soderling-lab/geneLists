#!/usr/bin/env Rscript

## Building a developmental brain disorder gene set
# collection for analysis with the AnRichment package.

# Data were scraped from: 
# "https://www.dbdb.urmc.rochester.edu/associations/list"
# To scrape the data, run ./scrape-UMRC-DBDB.py

## User parameters:
diseases <- c("autism",
	      "intellectual disability",
	      "attention deficit hyperactivity disorder",
	      "schizophrenia", 
	      "bipolar disorder",
	      "epilepsy")

# Imports.
suppressPackageStartupMessages({
	library(data.table)
	library(dplyr)
	library(getPPIs)
})

# Load functions in R/
devtools::load_all()

# Directories.
here <- getwd()
root <- dirname(dirname(here))
datadir <- file.path(root,"data")
rdatdir <- file.path(root,"rdata")
tabsdir <- file.path(root,"tables")
downdir <- file.path(root,"downloads")

# Load the raw data.
myfile <- file.path(downdir,"rochester-dbdb-associations.csv")
raw_data <- fread(myfile,drop=1)

# Map human gene symbols to Entrez.
genes <- unique(raw_data$Gene)
entrez <- mapIds(genes,from="symbol",to="entrez",species="human")
names(entrez) <- genes

# Map missing Entrez IDs by hand.
not_mapped <- genes[is.na(entrez)]
mapping_table <- fread(file.path(downdir,"not_mapped.txt"))
mapped_by_manual_search <- mapping_table$entrez
names(mapped_by_manual_search) <- mapping_table$gene
entrez[not_mapped] <- mapped_by_manual_search[names(entrez[not_mapped])]

# Check.
check <- sum(is.na(entrez)) == 0
if (!check)  { stop() }

# Add entrez IDs to data.
idy <- "Gene" # Column after which Entrez ids will be added.
data <- tibble::add_column(raw_data,
			   "hsEntrez" = entrez[raw_data[[idy]]],.after=idy)

# Map human genes in to their mouse homologs.
hsEntrez <- data$hsEntrez
msEntrez <- getHomologs(hsEntrez,taxid=10090)
data <- tibble::add_column(data,msEntrez=msEntrez,.after="hsEntrez")
# Remove rows with unmapped genes.
data <- data[msEntrez != "NA"] 

# Fix missing phenotype annotation for GNAQ and GNAS.
data$Phenotype[data$Gene=="GNAQ"] <- "Epilepsy;Intellectual disability"
data$Phenotype[data$Gene=="GNAS"] <- "Endocrine dysfunction;Intellectual disability"

# Split rows.
data <- tidyr::separate_rows(data,Phenotype,sep=";")

# Find genes associated with diseases of interest.
rows <- lapply(diseases,function(x) grep(x,tolower(data$Phenotype)))

# Split into disorder groups.
data_list <- lapply(rows,function(idx) data[idx,])

# Check disorder group sizes.
sizes <- sapply(data_list,function(x) length(unique(x$msEntrez)))

# Remove groups if size is 0.
data_list <- data_list[-which(sizes==0)]

# Status report.
nGenes <- length(unique(c(unlist((sapply(data_list,function(x) x$Gene))))))
nDisorders <- length(data_list)
message(paste("Compiled",nGenes,"mouse genes associated with",nDisorders,"DBDs!"))

# Write as gmt.
gmt_list <- lapply(data_list,function(x) x$msEntrez)
gmt_source = "https://www.dbdb.urmc.rochester.edu/associations/list"
gmt_file <- file.path(datadir,"5_URMC_DBDB.gmt")
write_gmt(gmt_list,gmt_source,gmt_file)
