#!/usr/bin/env Rscript

## Building a ASD-risk gene dataset from the WES data from 
# Sanders et al., 2015.

# Supplemental data were downloaded from:
# https://www.sciencedirect.com/science/article/pii/S0896627315007734?via%3Dihub#app2

## User parameters:
gene_source = "Sanders_et_al_2015"
data_source = "https://www.sciencedirect.com/science/article/pii/S0896627315007734?via%3Dihub#app2"
data_description = "ASD-associated genes."
output_file = "mouse_Sanders_ASD_geneSet.RData"

# Imports.
suppressPackageStartupMessages({
	library(anRichment)
	library(readxl)
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

# Load the Satterstrom data.
myfiles <- list.files(downdir,pattern="Sanders",full.names=TRUE)
names(myfiles) <- gsub(".xlsx","",sapply(strsplit(myfiles,"_"),"[",5))

# We need the data from S2. 
myfile <- myfiles["S7"]
raw_data <- read_excel(myfile)

# Keep genes with FDR < 0.1 from these columns.
"65genes_tadaFdrAscSscExomeSscAgpSmallDel"
"59genes_tadaFdrAscSscExome"

# Get human entrez ids.
hsEntrez <- data$entrez_id

# Status.
n_genes <- sum(!is.na(hsEntrez))
check <- n_genes == 102
if (!check) { stop("Warning, problem parsing excel data.") }

# Map human genes in to their mouse homologs.
msEntrez <- getHomologs(hsEntrez,taxid=10090)
data <- tibble::add_column(data,msEntrez=msEntrez,.after="entrez_id")

# Status.
n_mapped <- sum(!is.na(data$msEntrez))
message(paste(n_mapped,"of", n_genes, "human ASD-risk",
	      "genes were successfully mapped to mouse homologs."))

# Remove rows with unmapped genes.
data <- subset(data, !is.na(msEntrez))

# Split rows with multiple OMIM phenotypes.
data <- data %>% tidyr::separate_rows(OMIM_Phenotype,sep=";")

# Fill in blanks.
is_blank <- data$OMIM_Phenotype == "."
data$OMIM_Phenotype[is_blank] <- data$ASD_vs_DDID[is_blank]

# Save data.
fwrite(data,file.path(tabsdir,"Satterstrom_et_al_ASD_Genes.csv"))

# Split into disorder groups.
disorders <- unique(data$ASD_vs_DDID)
data_list <- data %>% group_by(ASD_vs_DDID) %>% group_split()
names(data_list) <- disorders

# Add combined data as well.
#data_list[["ASD-DDID"]] <- data

# Build gene sets:
geneSets <- list()
for (i in seq_along(data_list)) {
	subdat <- data_list[[i]]
	id <- paste0("Satterstrom_",names(data_list)[i])
	geneSets[[i]] <- newGeneSet(geneEntrez = unique(subdat$msEntrez),
				    geneEvidence = "IEA",
				    geneSource = gene_source, 
				    ID = id,
				    name = subdat$OMIM_Phenotype, # disorder name
				    description = data_description,
				    source = data_source,
				    organism = "mouse",
				    internalClassification = c("PL","Satterstrom_ASD"),
				    groups = "PL",
				    lastModified = Sys.Date())
} # Ends loop.

# Define gene collection groups.
PLgroup <- newGroup(name = "PL", 
		   description = data_description,
		   source = data_source)

# Combine as gene collection.
GOcollection <- newCollection(dataSets = geneSets, groups = list(PLgroup))

# Save as RData.
myfile <- file.path(rdatdir,output_file)
saveRDS(GOcollection,myfile)
