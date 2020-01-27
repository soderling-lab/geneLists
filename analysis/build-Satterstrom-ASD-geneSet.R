#!/usr/bin/env Rscript

## Building a ASD-risk gene dataset from the WES data from Satterstrom et al.,
## 2020.

# Supplemental data were downloaded from:
# https://www.cell.com/cell/fulltext/S0092-8674(19)31398-4?rss=yes#secsectitle0385

## User parameters:
gene_source = "Satterstrom_et_al_2020"
data_source = "https://www.biorxiv.org/content/10.1101/484113v3"
data_description = "ASD/DDID-associated genes."
output_file = "mouse_Satterstrom_ASD_geneSet.RData"

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
myfiles <- list.files(downdir,pattern="Satterstrom",full.names=TRUE)
names(myfiles) <- paste0("S",c(1:length(myfiles)))

# We need the data from S2. 
raw_data <- list()
sheet_names <- excel_sheets(myfiles["S2"])
for (sheet in sheet_names) {
	raw_data[[sheet]] <- read_excel(myfiles["S2"],sheet)
}

# Get ASD genes.
data <- raw_data[["102_ASD"]]

# Remove the last three rows--these just contain some comments.
data <- data[!is.na(data$entrez_id),]

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
#data <- data %>% tidyr::separate_rows(OMIM_Phenotype,sep=";")

# Fill in blanks.
#is_blank <- data$OMIM_Phenotype == "."
#data$OMIM_Phenotype[is_blank] <- data$ASD_vs_DDID[is_blank]

# Save data.
fwrite(data,file.path(tabsdir,"Satterstrom_et_al_ASD_Genes.csv"))

# Split into disorder groups.
#disorders <- unique(data$ASD_vs_DDID)
#data_list <- data %>% group_by(ASD_vs_DDID) %>% group_split()
#names(data_list) <- disorders

# Build gene sets:
geneSets <- newGeneSet(geneEntrez = unique(data$msEntrez),
		       geneEvidence = "IEA",
		       geneSource = gene_source, 
		       ID = "Satterstrom-ASD",
		       name = "Satterstrom-ASD",
		       description = data_description,
		       source = data_source,
		       organism = "mouse",
		       internalClassification = c("PL","Satterstrom_ASD"),
		       groups = "PL",
		       lastModified = Sys.Date())

# Define gene collection groups.
PLgroup <- newGroup(name = "PL", 
		   description = data_description,
		   source = data_source)

# Combine as gene collection.
DBDcollection <- newCollection(dataSets = list(geneSets), groups = list(PLgroup))

# Save as RData.
myfile <- file.path(rdatdir,output_file)
saveRDS(DBDcollection,myfile)
