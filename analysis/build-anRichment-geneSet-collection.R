#!/usr/bin/env Rscript

## Building a user-defined gene set for use with the anRichment package.
#  Adapted from the tutorial online:
#  https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/Tutorials/anRichment-Tutorial1.pdf

suppressPackageStartupMessages({
library(anRichment)
library(data.table)
library(dplyr)
})

# Download DisgeneNet currated gene list.
url <- "https://www.disgenet.org/static/disgenet_ap1/files/downloads/curated_gene_disease_associations.tsv.gz"
download.file(url,destfile=basename(url))

# Unzip.
system(command = paste("gunzip",basename(url)))

# Read with fread.
myfile <- tools::file_path_sans_ext(basename(url))
data <- fread(myfile)

# Remove raw data file.
unlink(myfile)

# Check out current evidence codes.
#knownEvidenceCodes()

# Split data into groups.
data_list <- split(data,data$diseaseId)

# Remove disease groups with less than 5 genes.
out <- seq_along(data_list)[sapply(data_list,function(x) dim(x)[1]<5)]
data_list <- data_list[-out]

# Loop to build gene sets:
geneSets <- list()

for (i in 1:length(data_list)) {
	df <- data_list[[i]]
	id <- names(data_list)[i]
	namen <- gsub(" ","_",df$diseaseName)
	geneSets[[i]] <- newGeneSet(geneEntrez = df$geneId,
				    geneEvidence = "IEA", # Inferred from Electronic Annotation
				    geneSource = "disgenet",
				    ID = id, # diseaseId
				    name = namen, # Shortened disease name
				    description = "Gene-disease associations from UNIPROT, CGI, ClinGen, Genomics England, CTD (human subset), PsyGeNET, and Orphanet.",
				    source = "https://www.disgenet.org/static/disgenet_ap1/files/downloads/curated_gene_disease_associations.tsv.gz",
				    organism = "human",
				    internalClassification = "DisGeneNET",
				    groups = "DisGeneNET",
				    lastModified = "2020-01-03")
}

# Define group.
PLgroup <- newGroup(name = "DisGeneNET", 
		   description = "Currated gene-disease associations from DisGenNET",
		   source = "disgenenet.org")
# Combine as gene collection.
DisGeneNETcollection <- newCollection(dataSets = geneSets, groups = list(PLgroup))

# Save.
myfile <- "DisGeneNETcollection.RData"
saveRDS(DisGeneNETcollection,myfile)
