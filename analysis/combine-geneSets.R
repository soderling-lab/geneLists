#!/usr/bin/env Rscript

## Building a developmental brain disorder gene set
# collection for analysis with the AnRichment package.

# Data download from: 
# "https://www.dbdb.urmc.rochester.edu/associations/list"

## User parameters:
map2mouse <- TRUE
save_data <- TRUE
filter_groups = TRUE
min_size = 3

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

# Load the raw data.
myfile <- file.path(downdir,"rochester-dbdb-associations.csv")
raw_data <- fread(myfile,drop=1)

# Map human gene symbols to Entrez.
genes <- unique(raw_data$Gene)
entrez <- mapIDs(genes,from="symbol",to="entrez",species="human")
names(entrez) <- genes

# Map missing Entrez IDs by hand.
message("Manually mapping missing ids...")
not_mapped <- genes[is.na(entrez)]
mapping_table <- fread(file.path(downdir,"not_mapped.txt"))
mapped_by_manual_search <- mapping_table$entrez
names(mapped_by_manual_search) <- mapping_table$gene
entrez[not_mapped] <- mapped_by_manual_search[names(entrez[not_mapped])]

# Check.
check <- sum(is.na(entrez)) == 0
if (check)  { message("Successfully mapped all human gene symbols to Entrez!") }

# Add entrez IDs to data.
idy <- "Gene" # Column after which Entrez ids will be added.
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

# Fix missing phenotype annotation for GNAQ and GNAS.
data$Phenotype[data$Gene=="GNAQ"] <- "Epilepsy;Intellectual disability"
data$Phenotype[data$Gene=="GNAS"] <- "Endocrine dysfunction;Intellectual disability"

# Split into disorder groups.
disorders <- unique(data$Phenotype)
data_list <- data %>% group_by(Phenotype) %>% group_split()
names(data_list) <- gsub(" ","_",disorders)

# Check disorder group sizes.
sizes <- sapply(data_list,function(x) length(unique(x$osEntrez)))

# Remove groups with less than min genes.
if (filter_groups) {
	keep <- seq(sizes)[sizes > min_size]
	nout <- length(data_list) - length(keep)
	data_list <- data_list[keep] # This limits to 15 disease groups.
	message(paste0("Number of disease groups removed with less than ",
		      min_size," genes: ",nout,"."))
}

# Status report.
nGenes <- length(unique(c(unlist((sapply(data_list,function(x) x$Gene))))))
nDisorders <- length(data_list)
message(paste("Compiled",nGenes,"genes associated with",nDisorders,"DBDs!"))

# Build gene sets:
geneSets <- list()
for (i in seq_along(data_list)) {
	subdat <- data_list[[i]]
	id <- paste0("DBDB-",names(data_list)[i])
	geneSets[[i]] <- newGeneSet(geneEntrez = unique(subdat$osEntrez),
				    geneEvidence = "IEA",
				    geneSource = "URMC-DBDB",
				    ID = id,
				    name = subdat$Phenotype[1], # disorder name
				    description = "URMC-DBD-associated genes.",
				    source = "https://www.dbdb.urmc.rochester.edu/associations/list",
				    organism = org,
				    internalClassification = c("PL","URMC-DBDB"),
				    groups = "PL",
				    lastModified = Sys.Date())
} # Ends loop.

# Define gene collection groups.
PLgroup <- newGroup(name = "PL", 
		   description = "G2P-DBD-associated genes.",
		   source = "https://www.ebi.ac.uk/gene2phenotype/")

# Combine as gene collection.
GOcollection <- newCollection(dataSets = geneSets, groups = list(PLgroup))

# Save as RData.
myfile <- file.path(rdatdir,paste0(org,"_","Rochester_URMC_DBDB",".RData"))
saveRDS(GOcollection,myfile)
