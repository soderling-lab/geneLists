#!/usr/bin/env Rscript

## Building a gene set 
# collection for analysis with the AnRichment package.

# Data download from: 
# https://www.ebi.ac.uk/gene2phenotype/downloads

## User parameters:
dataset <- "DD" # One of c("Cancer","DD","Eye","Skin")
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

# Functions.
devtools::load_all()

# Directories.
here <- getwd()
root <- dirname(dirname(here))
datadir <- file.path(root,"data")
rdatdir <- file.path(root,"rdata")
tabsdir <- file.path(root,"tables")
downdir <- file.path(root,"downloads")

# Downoad the raw data.
base_url <- "https://www.ebi.ac.uk/gene2phenotype/downloads"
datasets <- c(Cancer = "CancerG2P",
	      DD = "DDG2P",
	      Eye = "EyeG2P",
	      Skin = "SkinG2P")
myurl <- file.path(base_url,paste0(datasets[dataset],".csv.gz"))
download.file(myurl,destfile=basename(myurl),quiet=FALSE)
system(command = paste("gunzip",basename(myurl)))
myfile <- paste0(datasets[dataset],".csv")
raw_data <- fread(myfile)
unlink(myfile)

# Fix column names.
colnames(raw_data) <- gsub(" ", "_", colnames(raw_data))

# Map human gene symbols to Entrez.
genes <- unique(raw_data$gene_symbol)
nHsGenes <- length(unique(genes))
entrez <- mapIds(genes,from="symbol",to="entrez",species="human")
names(entrez) <- genes

# Map missing Entrez IDs by hand.
not_mapped <- genes[is.na(entrez)]
mapped_by_manual_search <- c('KIF1BP' = 26128,
			        'KARS' = 3735,
		                'MT-TP' = 4571,
		            	'QARS' = 5859,
				'HARS' = 3035,
				'HIST3H3' = 8290,
				'AARS' = 16,
				'ARSE' = 415,
				'DARS' = 1615,
				'HIST1H4B' = 8366,
				'IARS' = 3376,
				'HIST1H4C' = 8364,
				'HIST1H1E' = 3006,
				'HIST1H4J' = 8363,
				'EPRS' = 2058,
				'CARS' = 833,
				'TARS' = 689,
				'HIST1H2AC' = 8334)
entrez[not_mapped] <- mapped_by_manual_search[names(entrez[not_mapped])]

# Check.
check <- sum(is.na(entrez)) == 0
if (!check)  { stop() }

# Add entrez IDs to data.
idy <- "gene_symbol" # Column after which Entrez ids will be added.
data <- tibble::add_column(raw_data,"Entrez" = entrez[raw_data[[idy]]],.after=idy)

# Map human genes in to their mouse homologs.
hsEntrez <- data$Entrez
msEntrez <- getHomologs(hsEntrez,taxid=10090)
data <- tibble::add_column(data,msEntrez=msEntrez,.after="Entrez")
# Remove rows with unmapped genes.
data <- data[msEntrez != "NA"] 

# Get disorders that affect the brain/cognition.
disease_types <- "Brain/Cognition"
data <- tidyr::separate_rows(data,organ_specificity_list,sep=";")
data <- subset(data,data$organ_specificity_list==disease_types)

# Collect diseases of interest.
rows <- lapply(diseases,function(x) grep(x,tolower(data$disease_name)))

# Split into groups.
data_list <- lapply(rows,function(idx) data[idx,])

# Check disorder group sizes.
sizes <- sapply(data_list,function(x) length(unique(x$msEntrez)))

# Remove groups with 0 genes.
data_list <- data_list[-which(sizes==0)]

# Status report.
nGenes <- length(unique(unlist(sapply(data_list,function(x) x$msEntrez))))
nDisorders <- length(data_list)
message(paste("Compiled",nGenes,"mouse genes associated with",nDisorders,"DBDs!"))

# save as gmt.
gmt_list <- lapply(data_list,function(x) x$msEntrez)
gmt_source <- "https://www.ebi.ac.uk/gene2phenotype/downloads"
gmt_file <- file.path(datadir,"6_gene2phenotype.gmt")
write_gmt(gmt_list,gmt_source,gmt_file)
