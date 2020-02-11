#!/usr/bin/env Rscript

## Building a disease-associated gene GO collection for analysis with 
# the AnRichment package.

# Which dataset?
dataset <- "SFARI" # SFARI or Animal

# Imports.
suppressPackageStartupMessages({
	library(data.table)
	library(dplyr)
	library(getPPIs)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
datadir <- file.path(root,"data")
rdatdir <- file.path(root,"rdata")
tabsdir <- file.path(root,"tables")
downdir <- file.path(root,"downloads")

# Load functions.
devtools::load_all()

# Load the data.
datasets <- c(SFARI = "SFARI-Gene_genes",
	      Animal = "SFARI-Gene_animal-genes")
myfile <- list.files(downdir,pattern=datasets[dataset],full.names=TRUE)
data <- data.table::fread(myfile)
colnames(data) <- gsub("-","_",colnames(data))

# Map human gene symbols to Entrez.
genes <- data$gene_symbol
entrez <- mapIds(genes,from="symbol",to="entrez",species="human")
data <- tibble::add_column(data,"entrez_id" = entrez, .after = 4)
data <- data[entrez_id != "NA"]

# Map human genes in to their mouse homologs.
hsEntrez <- data$entrez_id
nHsGenes <- length(unique(hsEntrez))
msEntrez <- getHomologs(hsEntrez,taxid=10090)
data <- tibble::add_column(data,msEntrez=msEntrez,.after=5)
# Remove rows with unmapped genes.
data <- data[msEntrez != "NA"] 

# Status report.
nGenes <- length(unique(data$msEntrez))
message(paste0("Compiled ",nGenes," mouse genes associated with ",
	       "Autism spectrum disorders!"))

# Write data to file.
myfile <- file.path(tabsdir,paste0("mouse_",datasets[dataset],".csv"))
data.table::fwrite(data,myfile)

# Write to gmt.
gmt_list <- data$msEntrez
gmt_file <- file.path(datadir,"3_mouse_SFARI_Gene.gmt")
gmt_source <- "https://gene.sfari.org/tools/"
write_gmt(gmt_list,gmt_source,gmt_file)
