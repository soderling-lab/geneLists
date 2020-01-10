#!/usr/bin/env Rscript

## Building a gene set 
# collection for analysis with the AnRichment package.

# Data download from: 
# https://www.ebi.ac.uk/gene2phenotype/downloads

## User parameters:
dataset <- "DD" # One of c("Cancer","DD","Eye","Skin")
map2mouse <- TRUE
save_data <- TRUE
filter_groups <- FALSE
min_size <- 3

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

# Which organism, human or mouse?
if (map2mouse) { org <- "mouse" } else { org <- "human" }

# Downoad the raw data.
message("Downloading gene to phenotype database... \n")
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
entrez <- mapIDs(genes,from="symbol",to="entrez",species="human")
names(entrez) <- genes

# Map missing Entrez IDs by hand.
message("Manually mapping missing ids...")
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
if (check)  { message("Successfully mapped all human gene symbols to Entrez!") }

# Add entrez IDs to data.
idy <- "gene_symbol" # Column after which Entrez ids will be added.
data <- tibble::add_column(raw_data,"Entrez" = entrez[raw_data[[idy]]],.after=idy)

# Map human genes in to their mouse homologs.
if (map2mouse) {
	message("Mapping human genes to their mouse homologs...\n")
	hsEntrez <- data$Entrez
	msEntrez <- getHomologs(hsEntrez,taxid=10090)
	data <- tibble::add_column(data,osEntrez=msEntrez,.after="Entrez")
	# Remove rows with unmapped genes.
	n_out <- sum(is.na(msEntrez))
	percent_removed <- round(100*(n_out/length(msEntrez)),2)
	message(paste0("Genes without mouse homology: ",n_out,
		      " (",percent_removed," %)."))
	data <- data[osEntrez != "NA"] 
} else {
	hsEntrez <- data$entrez_id
	data <- tibble::add_column(data,osEntrez=osEntrez,.after="Entrez")
}

# Get disorders that affect the brain/cognition.
data <- tidyr::separate_rows(data,organ_specificity_list,sep=";")
data <- subset(data,data$organ_specificity_list=="Brain/Cognition")

# Status report.
nGenes <- length(unique(data$osEntrez))
nDisorders <- length(unique(data$disease_name))
message(paste("Compiled",nGenes,"genes associated with",nDisorders,"DBDs!"))

# Write data to file.
if (save_data) {
	myfile <- file.path(tabsdir,paste0(org,"_",datasets[dataset]))
	data.table::fwrite(data,myfile)
}

# Split into disorder groups.
disorders <- unique(data$disease_name)
data_list <- data %>% group_by(disease_name) %>% group_split()
names(data_list) <- disorders

# Check disorder group sizes.
sizes <- sapply(data_list,function(x) length(unique(x$osEntrez)))

# Remove groups with less than min genes.
if (filter_groups) {
	keep <- seq(sizes)[sizes > min_size]
	data_list <- data_list[keep] # This limits to 15 disease groups.
}

# Status report.
nGenes <- length(unique(c(unlist((sapply(data_list,function(x) x$gene_symbol))))))
nDisorders <- length(data_list)
message(paste("Compiled",nGenes,"genes associated with",nDisorders,"DBDs!"))

# Build gene sets:
geneSets <- list()
for (i in seq_along(data_list)) {
	subdat <- data_list[[i]]
	id <- paste0("DBD-",names(data_list)[i])
	geneSets[[i]] <- newGeneSet(geneEntrez = unique(subdat$osEntrez),
				    geneEvidence = "IEA",
				    geneSource = datasets[dataset],
				    ID = id,
				    name = subdat$disease_name[i], # disorder name
				    description = "G2P-DBD-associated genes.",
				    source = "https://www.ebi.ac.uk/gene2phenotype/downloads",
				    organism = org,
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
myfile <- file.path(rdatdir,paste0(org,"_",datasets[dataset],"_DBD",".RData"))
saveRDS(DBDcollection,myfile)
