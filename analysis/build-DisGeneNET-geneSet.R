#!/usr/bin/env Rscript

## Building a disease-associated gene GO collection for analysis with 
# the AnRichment package.

# Which DisGene dataset to download and compile.
# Choose one of: c("Curated_Disease_Genes","All_Disease_Genes",
# "Curated_Variants", "All_Variants")
#dataset <- "Curated_Variants"
dataset <- "All_Disease_Genes"
save_data <- TRUE # Save data to .csv?
map2mouse <- TRUE # Map human genes to mouse? If FALSE, then human data will be saved.

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

# Which organism?
if (map2mouse) {org <- "mouse" } else { org = "human" }

# Download and unzip the DisgeneNet currated gene list.
message("Downloading Disease-Gene associations from DisGeneNet.org...\n")
base_url <- "https://www.disgenet.org/static/disgenet_ap1/files/downloads"
datasets <- c(Curated_Disease_Genes = "curated_gene_disease_associations.tsv.gz",
	      All_Disease_Genes = "all_gene_disease_associations.tsv.gz",
	      Curated_Variants = "curated_variant_disease_associations.tsv.gz",
	      All_Variants = "all_variant_disease_associations.tsv.gz",
	      Variant_Gene_Map = "variant_to_gene_mappings.tsv.gz")
myurl <- file.path(base_url,datasets[dataset])
download.file(myurl,destfile=basename(myurl), quiet=FALSE)
system(command = paste("gunzip",basename(myurl)))
# Read the data with fread.
myfile <- tools::file_path_sans_ext(basename(myurl))
data <- fread(myfile)
# Remove raw data file.
unlink(myfile)

# If working with gene variants, get mapping data.
if (grepl("Variants",dataset)) {
	# Download.
	message("Downloading Variant-Gene mapping data from DisGeneNet.org...\n")
	myurl <- file.path(base_url,datasets["Variant_Gene_Map"])
	download.file(myurl,destfile=basename(myurl),quiet=FALSE)
	system(command = paste("gunzip",basename(myurl)))
	# Read the data with fread.
	myfile <- tools::file_path_sans_ext(basename(myurl))
	varmap <- fread(myfile)
	# Add geneID and geneSymbol columns to data.
	entrez <- varmap$geneId[match(data$snpId,varmap$snpId)]
	genes <- varmap$geneSymbol[match(data$snpId,varmap$snpId)]
	data <- tibble::add_column(data,geneId=entrez,.after=1)
	data <- tibble::add_column(data,geneSymbol=genes,.after=2)
	# Remove unmapped variants.
	data <- data[geneId != "NA"]
	# Remove raw data file.
	unlink(myfile)
}

# DisGeneNET is a database of human genes and their associated diseases.
# Map human genes in to their mouse homologs.
if (map2mouse) {
	message("Mapping human genes to their mouse homologs...\n")
	hsEntrez <- data$geneId
	msEntrez <- getHomologs(hsEntrez,taxid=10090)
	data <- tibble::add_column(data,osEntrez=msEntrez,.after=1)
	# Remove rows with unmapped genes.
	n_out <- sum(is.na(msEntrez))
	percent_removed <- round(100*(n_out/length(msEntrez)),2)
	message(paste("Percent disease genes without mouse homology:",
		      percent_removed))
	data <- subset(data, !is.na(data$osEntrez))
} else {
	hsEntrez <- data$geneId
	data <- tibble::add_column(data,osEntrez=osEntrez,.after=1)
}

# Write data to file.
if (save_data) {
	myfile <- file.path(tabsdir,paste0(org,"_DisGeneNet_",dataset,".csv"))
	data.table::fwrite(data,myfile)
}

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
message("Building disease gene sets...\n")
geneSets <- list()
pbar <- txtProgressBar(min=1,max=length(data_list),style=3)
for (i in 1:length(data_list)) {
	setTxtProgressBar(pbar,i)
	df <- data_list[[i]]
	diseaseID <- names(data_list)[i]
	entrez <- df$osEntrez
	namen <- gsub(" ","_",unique(df$diseaseName))
	if (length(namen) > 1) { stop() }
	geneSets[[i]] <- newGeneSet(geneEntrez = entrez,
				    geneEvidence = "IEA", # Inferred from Electronic Annotation
				    geneSource = "DisGeneNET",
				    ID = diseaseID,
				    name = namen, # disease name
				    description = data_description,
				    source = myurl,
				    organism = org,
				    internalClassification = c("PL","DisGeneNET"),
				    groups = "PL",
				    lastModified = "2020-01-03")
	if (i==length(data_list)) { close(pbar); message("\n") }
} # Ends loop.

# Define gene collection groups.
PLgroup <- newGroup(name = "PL", 
		   description = "Gene-disease associations from DisGenNET.",
		   source = "disgenenet.org")

# Combine as gene collection.
DisGeneNETcollection <- newCollection(dataSets = geneSets, groups = list(PLgroup))

# Save as RData.
myfile <- file.path(rdatdir,paste0(org,"_DisGeneNet_",dataset,".RData"))
saveRDS(DisGeneNETcollection,myfile)
