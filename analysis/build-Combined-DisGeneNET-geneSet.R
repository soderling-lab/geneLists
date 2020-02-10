#!/usr/bin/env Rscript

## Building a disease-associated gene GO collection for analysis with 
# the AnRichment package.

# Which DisGene dataset to download and compile.
dataset <- "All_Variants"
disease_types <- c("Mental or Behavioral Dysfunction", "Mental Process")
min_size <- 3
max_size <- 500 
# Dont include BEFREE and CTD databases.
data_sources <- c("PSYGENET","LHGDN","HPO",
		      "GENOMICS_ENGLAND","RGD","GWASCAT","GWASDB","MGD",
		      "CLINVAR","UNIPROT","ORPHANET")

# All potential data sources:
all_data_sources <- c("CTD_human","BEFREE","PSYGENET","LHGDN","HPO",
		      "GENOMICS_ENGLAND","RGD","GWASCAT","GWASDB","MGD",
		      "CLINVAR","UNIPROT","CTD_rat","ORPHANET","CTD_mouse")
## Datasets:
# [1] Curated_Disease_Genes
# [2] All_Disease_Genes
# [3] Curated_Variants
# [4] All_Variants

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

# Notin function.
`%notin%` <- Negate(`%in%`)

# Url for downloads.
base_url <- "https://www.disgenet.org/static/disgenet_ap1/files/downloads"

# Define a function that gets DisGeneNet data:
getDisGeneNetData <- function(dataset,disease_types,data_sources,base_url){
	# Download, unzip, and load the DisGeneNet data.
	datasets <- c(Curated_Disease_Genes = "curated_gene_disease_associations.tsv.gz",
		      All_Disease_Genes = "all_gene_disease_associations.tsv.gz",
		      Curated_Variants = "curated_variant_disease_associations.tsv.gz",
		      All_Variants = "all_variant_disease_associations.tsv.gz",
		      Variant_Gene_Map = "variant_to_gene_mappings.tsv.gz")
	myurl <- file.path(base_url,datasets[dataset])
	download.file(myurl,destfile=basename(myurl), quiet=TRUE)
	system(command = paste("gunzip",basename(myurl)))
	myfile <- tools::file_path_sans_ext(basename(myurl))
	data <- fread(myfile)
	unlink(myfile) # Remove raw data file.
	# If working with gene variants, get mapping data from DisGeneNet.
	if (grepl("Variants",dataset)) {
		myurl <- file.path(base_url,datasets["Variant_Gene_Map"])
		download.file(myurl,destfile=basename(myurl),quiet=TRUE)
		system(command = paste("gunzip",basename(myurl)))
		myfile <- tools::file_path_sans_ext(basename(myurl))
		varmap <- fread(myfile)
		# Add geneID and geneSymbol columns to data.
		entrez <- varmap$geneId[match(data$snpId,varmap$snpId)]
		genes <- varmap$geneSymbol[match(data$snpId,varmap$snpId)]
		data <- tibble::add_column(data,geneId=entrez,.after=1)
		data <- tibble::add_column(data,geneSymbol=genes,.after=2)
		data <- data[geneId != "NA"] # Remove unmapped variants.
		unlink(myfile) # Remove raw data file.
	}
	# DisGeneNET is a database of human genes and their associated diseases.
	# Map human genes in to their mouse homologs.
	hsEntrez <- data$geneId
	msEntrez <- getHomologs(hsEntrez,taxid=10090)
	nMsGenes <- length(unique(msEntrez))
	data <- tibble::add_column(data,msEntrez=msEntrez,.after=1)
	# Remove rows with unmapped genes.
	data <- subset(data, !is.na(data$msEntrez))
	# Get diseases realated to brain function.
	data <- subset(data,data$diseaseSemanticType %in% disease_types)
	# Split rows (genes) derived from multiple sources.
	data <- data %>% tidyr::separate_rows(source,sep=";")
	# Keep the data of interest.
	data <- subset(data,data$source %in% data_sources)
	return(data)
}

# Get all DisGeneNet datasets.
datasets <- c("Curated_Disease_Genes","All_Disease_Genes","Curated_Variants","All_Variants")
disease_types <- c("Mental or Behavioral Dysfunction", "Mental Process")
all_data <- lapply(datasets,function(x) getDisGeneNetData(x,disease_types,data_sources,base_url))
names(all_data) <- datasets
backup_data <- all_data

## Combine datasets.
colNames <- intersect(colnames(all_data[[1]]),colnames(all_data[[4]]))
allNames <- union(colnames(all_data[[1]]),colnames(all_data[[4]]))

# Note: The following columns will be dropped:
#allNames[allNames %notin% colNames]

# Merge dataframes.
data <- all_data %>% purrr::reduce(full_join,by=colNames)

# If data was found in multiple dataframes a column was added. 
# e.g. chromosome.x; remove these columns.
data <- data[,-c(grep("\\.",colnames(data)))]

# Remove non-unique rows.
data <- unique(data)

# Split data into disease groups.
data_list <- split(data,data$diseaseId)

# Check disorder group sizes.
sizes <- sapply(data_list,function(x) length(unique(x$msEntrez)))

# Remove groups with less than min genes or greater than max genes.
keep <- names(sizes)[sizes > min_size & sizes < max_size]
data_list <- data_list[keep] 
data <- do.call(rbind,data_list)

# Status report.
nGenes <- length(unique(data$msEntrez))
nDisorders <- length(unique(data$diseaseId))
message(paste0("Compiled ",nGenes," mouse genes associated with ",
	       nDisorders," DBDs or phenotypes!"))

# Save to file.
myfile <- file.path(tabsdir,paste0("mouse_Combined_DisGeneNet_geneSet.csv"))
fwrite(data,myfile)

# Description of the data.
data_description <- paste("Gene-disease associations compiled from",
			  paste(data_sources,collapse=", "))

# Loop to build gene sets:
geneSets <- list()
for (i in 1:length(data_list)) {
	df <- data_list[[i]]
	diseaseID <- names(data_list)[i]
	entrez <- df$msEntrez
	namen <- gsub(" ","_",unique(df$diseaseName))
	geneSets[[i]] <- newGeneSet(geneEntrez = entrez,
				    geneEvidence = "IEA", # Inferred from Electronic Annotation
				    geneSource = paste("DisGeneNET",dataset,sep="-"),
				    ID = diseaseID,
				    name = namen, # disease name
				    description = data_description,
				    source = base_url,
				    organism = "mouse",
				    internalClassification = c("PL","DisGeneNET"),
				    groups = "PL",
				    lastModified = "2020-01-03")
} # Ends loop.

# Define gene collection groups.
PLgroup <- newGroup(name = "PL", 
		   description = "Gene-disease associations from DisGenNET.",
		   source = "disgenenet.org")

# Combine as gene collection.
DisGeneNETcollection <- newCollection(dataSets = geneSets, groups = list(PLgroup))

# Save as RData.
myfile <- file.path(rdatdir,paste0("mouse_Combined_DisGeneNet_geneSet.RData"))
saveRDS(DisGeneNETcollection,myfile)
