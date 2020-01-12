#!/usr/bin/env Rscript

## Building Synaptic GO collection for GO analysis with the 
# AnRichment package.

##  SynGO data downloaded from: https://www.syngoportal.org/ 

## User parameters:
release <- "20180731" # Most recent SynGO release.
map2mouse <- TRUE

# Imports.
suppressPackageStartupMessages({
	library(readxl)
	library(tidyr)
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
message("Downloading Synaptic GO (SyngGO) database... \n")
base_url <- "https://syngoportal.org/data/download.php?file=SynGO_bulk_download_release_"
myurl <- paste0(base_url,release,".zip")
myfile <- file.path(downdir,sapply(strsplit(basename(myurl),"file="),"[",2))
download.file(myurl,destfile=myfile,method="curl",quiet=TRUE,extra="--insecure")
cmd <- paste0("unzip -oq ",myfile," -d /",tools::file_path_sans_ext(myfile))
system(command = cmd)
unlink(myfile)

# Load SynGO data.
mydir <- tools::file_path_sans_ext(myfile)
myfiles <- list.files(mydir,pattern="xlsx",full.names=TRUE)
syngo_data <- sapply(myfiles,read_excel)
names(syngo_data) <- gsub(".xlsx","",basename(names(syngo_data)))  

# Get mouse Entrez IDs by mapping MGI IDs to Entrez.
genes_df <- syngo_data$syngo_genes 
# Seperate rows with multiple MGI ids.
# Warning caused by missing (NA) MGI ids can be ignored.
genes_df <- genes_df %>% separate_rows(mgi_id,sep=",")
mgi <- paste0("MGI:",genes_df$mgi_id)
msEntrez <- mapIDs(mgi,from="mgi",to="entrez",species="mouse")
names(msEntrez) <- genes_df$mgi_id
genes_df$msEntrez <- msEntrez
if (map2mouse){
	# Use mouse entrez.
	genes_df$osEntrez <- msEntrez
} else {
	# Use human entrez.
	genes_df$osEntrez <- genes_df$entrez_id
}

# Collect GO annotations and add this column to syngo genes.
anno_df <- syngo_data$syngo_annotations
idx <- match(anno_df$"human ortholog gene hgnc_id",genes_df$hgnc_id)
go_df <- data.frame("id" = as.character(anno_df$"GO term ID"), 
		    "domain" = as.character(anno_df$"GO domain"),
		    "pmid" = as.character(anno_df$PMID), 
		    "name" = as.character(anno_df$"GO term name"),
		    "osEntrez" = as.character(genes_df$osEntrez[idx]))

# Group into gene groups.
data_list <- go_df %>% group_by(id) %>% group_split()
names(data_list) <- as.character(sapply(data_list,function(x) unique(x$id)))

# Loop to build gene sets:
geneSets <- list()

for (i in seq_along(data_list)) {
	id <- names(data_list)[i]
	subdat <- data_list[[i]]
	geneSets[[i]] <- newGeneSet(geneEntrez = subdat$osEntrez,
				    geneEvidence = "IEA", 
				    geneSource = "SynGO",
				    ID = id, # diseaseId
				    name = unique(subdat$name), 
				    description = "Synaptic gene ontology",
				    source = myurl,
				    organism = org,
				    internalClassification = c("PL","syngo"),
				    groups = "PL",
				    lastModified = Sys.Date())
}

# Define PLgroup.
description <- "Currated synaptic gene ontology from SynGO database."
PLgroup = newGroup(name = "PL", 
		   description = description,
		   source = "syngoportal.org")

# Combine as gene collection.
SynGOcollection <- newCollection(dataSets = geneSets, groups = list(PLgroup))

# Save.
myfile <- file.path(rdatdir,paste0(org,"_SynGO_geneSet.RData"))
saveRDS(SynGOcollection,"mouse_SynGO_geneSet.RData")
