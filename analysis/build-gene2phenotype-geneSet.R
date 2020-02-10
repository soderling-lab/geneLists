#!/usr/bin/env Rscript

## Building a gene set 
# collection for analysis with the AnRichment package.

# Data download from: 
# https://www.ebi.ac.uk/gene2phenotype/downloads

## User parameters:
dataset <- "DD" # One of c("Cancer","DD","Eye","Skin")
min_size <- 3
max_size <- 500

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

# Split into disorder groups.
disorders <- unique(data$disease_name)
data_list <- data %>% group_by(disease_name) %>% group_split()
names(data_list) <- disorders

# Check disorder group sizes.
sizes <- sapply(data_list,function(x) length(unique(x$msEntrez)))

# Remove groups with less than min genes.
keep <- names(sizes)[ sizes > min_size & sizes < max_size]
data_list <- data_list[keep] # This limits to 15 disease groups.
data <- do.call(rbind,data_list)

# Status report.
nGenes <- length(unique(data$msEntrez))
nDisorders <- length(data_list)
message(paste("Compiled",nGenes,"mouse genes associated with",nDisorders,"DBDs!"))

# Save data to file.
myfile <- file.path(tabsdir,paste0("mouse_",datasets[dataset],"_DBD_geneSet.csv"))
fwrite(data,myfile)

# Build gene sets:
geneSets <- list()
for (i in seq_along(data_list)) {
	subdat <- data_list[[i]]
	id <- gsub(" ","_",names(data_list)[i])
	geneSets[[i]] <- newGeneSet(geneEntrez = unique(subdat$msEntrez),
				    geneEvidence = "IEA",
				    geneSource = datasets[dataset],
				    ID = id,
				    name = subdat$disease_name[1], 
				    description = "G2P-DBD-associated genes.",
				    source = "https://www.ebi.ac.uk/gene2phenotype/downloads",
				    organism = "mouse",
				    internalClassification = c("PL","G2P-DBD"),
				    groups = "PL",
				    lastModified = Sys.Date())
} # Ends loop.

# Define gene collection groups.
PLgroup <- newGroup(name = "PL", 
		   description = "G2P-DBD-associated genes.",
		   source = "https://www.ebi.ac.uk/gene2phenotype/")

# Combine as gene collection.
DBDcollection <- newCollection(dataSets = geneSets, groups = list(PLgroup))

# Save as RData.
myfile <- file.path(rdatdir,paste0("mouse_",datasets[dataset],"_DBD_geneSet.RData"))
saveRDS(DBDcollection,myfile)
