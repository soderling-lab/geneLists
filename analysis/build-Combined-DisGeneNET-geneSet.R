#!/usr/bin/env Rscript

## Building a disease-associated gene GO collection for analysis with 
# the AnRichment package.

## Parameters:
nThreads <- 7

## DisGeneNet Datasets:
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
	library(parallel)
	library(doParallel)
})

# Directories.
here <- getwd()
root <- dirname(here)
gmtdir <- file.path(root,"gmt")
rdatdir <- file.path(root,"rdata")
tabsdir <- file.path(root,"tables")

# Load functions.
devtools::load_all()  

# Misc function.
`%notin%` <- Negate(`%in%`)

# Define a function that gets DisGeneNet data:
getDisGeneNetData <- function(dataset = c("Curated_Disease_Genes","All_Disease_Genes",
					  "Curated_Variants","All_Variants"),
			      data_sources = c("CTD_human","BEFREE","PSYGENET","LHGDN","HPO",
					       "GENOMICS_ENGLAND","RGD","GWASCAT","GWASDB","MGD",
					       "CLINVAR","UNIPROT","CTD_rat","ORPHANET","CTD_mouse"),
			      base_url = "https://www.disgenet.org/static/disgenet_ap1/files/downloads") {
	# Download, unzip, and load the DisGeneNet data.
	# Variant Gene map is used for mapping gene variants.
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
	# Split rows (genes) derived from multiple sources.
	data <- data %>% tidyr::separate_rows(source,sep=";")
	# Keep the data of interest.
	data <- subset(data,data$source %in% data_sources)
	return(data)
}

# Get all DisGeneNet datasets.
datasets <- c("Curated_Disease_Genes","All_Disease_Genes","Curated_Variants","All_Variants")
all_data <- lapply(datasets,function(x) getDisGeneNetData(x))
names(all_data) <- datasets

## Combine datasets by shared column names.
colNames <- Reduce(intersect,sapply(all_data,colnames))

# Note: The following columns will be dropped:
#allNames <- Reduce(union,sapply(all_data,colnames))
#allNames[allNames %notin% colNames]

# Merge dataframes.
data <- all_data %>% purrr::reduce(full_join,by=colNames)

# If data was found in multiple dataframes then a column is added.
# e.g. chromosome.x; remove these duplicated columns.
data <- data[,-c(grep("\\.",colnames(data)))]

# Remove non-unique rows.
#data <- unique(data)

# We are interested in genes associated with:
diseases <- c("autism",
	      "intellectual disability",
	      "attention deficit hyperactivity disorder",
	      "schizophrenia", 
	      "bipolar disorder",
	      "epilepsy")

# Write a function to find these rows.
getRows <- function(diseases,i) {
	idx <- grep(diseases[i], tolower(data$diseaseName))
	return(idx)
}

# Parallelize the task with dopar.
workers <- makeCluster(nThreads)
registerDoParallel(workers)
rows <- foreach(i=seq_along(diseases)) %dopar% { getRows(diseases,i) }
suppressWarnings(stopCluster(workers))

# Split data into disease groups.
data_list <- lapply(rows, function(idx) data[idx,])
names(data_list) <- diseases

# Check disorder group sizes.
sizes <- sapply(data_list,function(x) length(unique(x$msEntrez)))

# Write as gmt file.
gmt_list <- lapply(data_list,function(x) x$msEntrez)
gmt_file <- file.path(gmtdir,"mouse_DisGeneNET.gmt")
write_gmt(gmt_list,gmt_source="DisGeneNET",gmt_file)

# Status report.
nGenes <- length(unique(unlist(sapply(data_list,function(x) x$msEntrez))))
nDisorders <- length(diseases)
message(paste0("Compiled ",nGenes," mouse genes associated with ",
	       nDisorders," DBDs!"))

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
