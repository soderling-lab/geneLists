#!/usr/bin/env Rscript

## Building a ASD-risk gene dataset from the WES data from 
# Sanders et al., 2015.

# Supplemental data were downloaded from:
# https://www.sciencedirect.com/science/article/pii/S0896627315007734?via%3Dihub#app2

## User parameters:
gene_source = "Sanders_et_al_2015"
data_source = "https://www.ncbi.nlm.nih.gov/pubmed/26402605"
data_description = "ASD-associated genes."
output_file = "mouse_Sanders_ASD_geneSet.RData"

# Imports.
suppressPackageStartupMessages({
	library(anRichment)
	library(readxl)
	library(data.table)
	library(dplyr)
	library(getPPIs)
})

# Directories.
here <- getwd()
root <- dirname(here)
funcdir <- file.path(root,"R")
rdatdir <- file.path(root,"rdata")
tabsdir <- file.path(root,"tables")
downdir <- file.path(root,"downloads")

# Functions.
invisible(sapply(list.files(funcdir,pattern=".R",full.names=TRUE),source))

# Load the data.
myfiles <- list.files(downdir,pattern="Sanders",full.names=TRUE)
names(myfiles) <- gsub(".xlsx","",sapply(strsplit(myfiles,"_"),"[",5))

# We need the data from S7.
myfile <- myfiles["S7"]
data <- read_excel(myfile)

# Keep genes with FDR < 0.1 from these columns.
cols <- c("65genes_tadaFdrAscSscExomeSscAgpSmallDel",
	  "59genes_tadaFdrAscSscExome")
keep <- data[cols[1]] != "0" | data[cols[2]] != "0"
data <- subset(data,keep)

# Map gene symbols to entrez.
symbols <- data$RefSeqGeneName
entrez <- mapIDs(symbols,from="symbol",to="entrez",species="human")

# Map missing by hand.
not_mapped <- which(is.na(entrez))
entrez[not_mapped] <- c(51111,55914)

# Add to data.
data <- tibble::add_column(data,"entrez"=entrez,.after="TadaGeneName")

# Map human genes in to their mouse homologs.
msEntrez <- getHomologs(entrez,taxid=10090)
data <- tibble::add_column(data,msEntrez=msEntrez,.after="entrez")

# Status.
n_mapped <- sum(!is.na(msEntrez))
n_genes <- length(entrez)
message(paste(n_mapped,"of", n_genes, "human ASD-risk",
	      "genes were successfully mapped to mouse homologs."))

# Save as gmt.
gmt_file <- file.path(rdatdir,"human_Sanders_ASD.gmt")
write_gmt(list("Sanders-ASD" = data$entrez),data_source,gmt_file)

# Remove rows with unmapped genes.
data <- subset(data, !is.na(msEntrez))

# Save as gmt.
gmt_file <- file.path(rdatdir,"mouse_Sanders_ASD.gmt")
write_gmt(list("Sanders-ASD" = data$msEntrez),data_source,gmt_file)

# Build gene set.
geneSets <- newGeneSet(geneEntrez = unique(data$msEntrez),
		       geneEvidence = "IEA",
		       geneSource = gene_source, 
		       ID = "Sanders-ASD",
		       name = "Sanders-ASD", # disorder name
		       description = data_description,
		       source = data_source,
		       organism = "mouse",
		       internalClassification = c("PL","Sanders_ASD"),
		       groups = "PL",
		       lastModified = Sys.Date())

# Define gene collection groups.
PLgroup <- newGroup(name = "PL", 
		   description = data_description,
		   source = data_source)

# Combine as gene collection.
DBDcollection <- newCollection(dataSets = list(geneSets), groups = list(PLgroup))

# Save as RData.
myfile <- file.path(rdatdir,output_file)
saveRDS(DBDcollection,myfile)
