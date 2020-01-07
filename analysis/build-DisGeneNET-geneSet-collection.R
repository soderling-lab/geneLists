#!/usr/bin/env Rscript

## Building a disease-associated gene GO collection for analysis with 
# the AnRichment package.

# Imports.
suppressPackageStartupMessages({
	library(anRichment)
	library(data.table)
	library(dplyr)
	library(getPPIs)
})

# Directories.
here <- getwd()
root <- dirname(here)
rdatdir <- file.path(root,"rdata")

# Download and unzip the DisgeneNet currated gene list.
myurl <- "https://www.disgenet.org/static/disgenet_ap1/files/downloads/curated_gene_disease_associations.tsv.gz"
download.file(myurl,destfile=basename(myurl))
system(command = paste("gunzip",basename(myurl)))

# Read the data with fread.
myfile <- tools::file_path_sans_ext(basename(myurl))
data <- fread(myfile)

# Remove raw data file.
unlink(myfile)

# DisGeneNET is a database of human genes and their associated diseases.
# Map human genes in to their mouse homologs.
hsEntrez <- data$geneId
msEntrez <- getHomologs(hsEntrez,taxid=10090)
data <- tibble::add_column(data,msEntrez,.after=1)

# Remove rows with unmapped genes.
n_out <- sum(is.na(msEntrez))
percent_removed <- round(100*(n_out/length(msEntrez)),2)
message(paste("Percent disease genes without mouse homology:",
	      percent_removed))
data <- subset(data, !is.na(data$msEntrez))

# Split data into disease groups.
data_list <- split(data,data$diseaseId)
nDiseases <- length(data_list)
message(paste("Number of disease groups:",nDiseases))

# Remove disease groups with less than 5 genes.
#out <- seq_along(data_list)[sapply(data_list,function(x) dim(x)[1]<5)]
#data_list <- data_list[-out]

# Description of the data.
data_description <- paste("Gene-disease associations from UNIPROT, CGI, 
			  ClinGen, Genomics England, CTD (human subset), 
			  PsyGeNET, and Orphanet.")

# Loop to build gene sets:
geneSets <- list()
for (i in 1:length(data_list)) {
	df <- data_list[[i]]
	diseaseID <- names(data_list)[i]
	entrez <- df$msEntrez
	namen <- gsub(" ","_",df$diseaseName)
	geneSets[[i]] <- newGeneSet(geneEntrez = entrez,
				    geneEvidence = "IEA", # Inferred from Electronic Annotation
				    geneSource = "DisGeneNET",
				    ID = diseaseID,
				    name = namen, # disease name
				    description = data_description,
				    source = myurl,
				    organism = "mouse",
				    internalClassification = "DisGeneNET",
				    groups = "DisGeneNET",
				    lastModified = "2020-01-03")
}

# Define gene collection groups.
PLgroup <- newGroup(name = "DisGeneNET", 
		   description = "Currated gene-disease associations from DisGenNET mapped to mouse homologs.",
		   source = "disgenenet.org")

# Combine as gene collection.
DisGeneNETcollection <- newCollection(dataSets = geneSets, groups = list(PLgroup))

# Save.
myfile <- file.path(rdatdir,"mouse_DisGeneNETcollection.RData")
saveRDS(DisGeneNETcollection,myfile)
