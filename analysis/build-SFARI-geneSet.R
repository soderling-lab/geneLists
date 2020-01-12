#!/usr/bin/env Rscript

## Building a disease-associated gene GO collection for analysis with 
# the AnRichment package.

# Which dataset?
dataset <- "SFARI" # SFARI or Animal
map2mouse <- TRUE
save_data <- TRUE

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
tabsdir <- file.path(root,"tables")
downdir <- file.path(root,"downloads")

# Which organism?
if (map2mouse) { org <- "mouse" } else { org <- "human" }

# Load the data.
datasets <- c(SFARI = "SFARI-Gene_genes",
	      Animal = "SFARI-Gene_animal-genes")
myfile <- list.files(downdir,pattern=datasets[dataset],full.names=TRUE)
data <- data.table::fread(myfile)
colnames(data) <- gsub("-","_",colnames(data))

# Map human gene symbols to Entrez.
genes <- data$gene_symbol
entrez <- mapIDs(genes,from="symbol",to="entrez",species="human")
data <- tibble::add_column(data,"entrez_id" = entrez, .after = 4)
data <- data[entrez_id != "NA"]

# Map human genes in to their mouse homologs.
if (map2mouse) {
	message("Mapping human genes to their mouse homologs...\n")
	hsEntrez <- data$entrez_id
	msEntrez <- getHomologs(hsEntrez,taxid=10090)
	data <- tibble::add_column(data,osEntrez=msEntrez,.after=5)
	# Remove rows with unmapped genes.
	n_out <- sum(is.na(msEntrez))
	percent_removed <- round(100*(n_out/length(msEntrez)),2)
	message(paste("Percent disease genes without mouse homology:",
		      percent_removed))
	data <- data[osEntrez != "NA"] 
} else {
	hsEntrez <- data$entrez_id
	data <- tibble::add_column(data,osEntrez=osEntrez,.after=5)
}

# Write data to file.
if (save_data) {
	myfile <- file.path(tabsdir,paste0(org,"_",datasets[dataset],".csv"))
	data.table::fwrite(data,myfile)
}

# Build gene set:
geneSet <- newGeneSet(geneEntrez = data$osEntrez,
		      geneEvidence = "IEA", # Inferred from Electronic Annotation
		      geneSource = datasets[dataset],
		      ID = datasets[dataset],
		      name = "SFARI", # disease name
		      description = "SFARI autism genes.",
		      source = "https://gene.sfari.org/tools/",
		      organism = org,
		      internalClassification = c("PL","SFARI"),
		      groups = "PL",
		      lastModified = Sys.Date())

# Define gene collection groups.
PLgroup <- newGroup(name = "PL", 
		   description = "SFARI autism-associated genes.",
		   source = "sfari.org")

# Combine as gene collection.
SFARIcollection <- newCollection(dataSets = list(geneSet), groups = list(PLgroup))

# Save as RData.
myfile <- file.path(rdatdir,paste0(org,"_",datasets[dataset],"_geneSet.RData"))
saveRDS(SFARIcollection,myfile)
