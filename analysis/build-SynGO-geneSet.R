#!/usr/bin/env Rscript

## Building Synaptic GO collection for GO analysis with the 
# AnRichment package.
#  SynGO downloaded from: https://www.syngoportal.org/ 

## User parameters:
release <- "20180731"
map2mouse <- TRUE
save_data <- TRUE

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
download.file(myurl,destfile=myfile,method="curl",quiet=FALSE,extra="--insecure")
cmd <- paste0("unzip ",myfile," -d /",tools::file_path_sans_ext(myfile))
system(command = cmd)
unlink(myfile)

# Load SynGO data.
mydir <- tools::file_path_sans_ext(myfile)
myfiles <- list.files(mydir,pattern="xlsx",full.names=TRUE)
data_list <- sapply(myfiles,read_excel)
names(data_list) <- gsub(".xlsx","",basename(names(data_list)))  


## Map mgi_ids to entrez.
# Seperate rows with multiple MGI ids.
# Warning caused by missing (NA) MGI ids can be ignored.
genes <- data_list$syngo_genes %>% separate_rows(mgi_id,sep=",")
mgi <- paste0("MGI:",genes$mgi_id)
entrez <- mapIDs(mgi,from="mgi",to="entrez",species="mouse")
names(entrez) <- genes$mgi_id
genes$mus_entrez <- entrez

# Collect GO annotations.
# Add column to syngo annotations.
anno_df <- data_list$syngo_annotations
idx <- match(anno_df$"human ortholog gene hgnc_id",genes$hgnc_id)
mus_entrez <- genes$mus_entrez[idx]
goDat <- data.frame("id" = paste0("SynGO:",anno_df$"SynGO ID"), 
		    "domain" = anno_df$"GO domain",
		    "pmid" = anno_df$PMID, 
		    "name" = anno_df$"GO term name",
		    "entrez" = mus_entrez)
goDat$entrez <- as.character(goDat$entrez)
goDat$domain <- as.character(goDat$domain)
goDat$name <- as.character(goDat$name)

# Split into two GO domains.
goDomains <- goDat %>% dplyr::group_by(domain) %>% dplyr::group_split()
names(goDomains) <- c("BP","CC")

# Collect BP and CC pathways into list.
bp <- goDomains[["BP"]]
cc <- goDomains[["CC"]]
pathways <- list("BP" = split(bp$entrez,bp$name),
		 "CC" = split(cc$entrez,cc$name))

# Save.
myfile <- file.path(rdatdir,"SynGO.RData")
saveRDS(pathways,myfile)

# Arbitrary names.
data_list <- do.call(c,pathways)
namen <- gsub(" ","_",names(data_list))
names(data_list) <- paste0("SYNGO_",c(1:length(data_list)))

# Loop to build gene setst:
geneSets <- list()
for (i in 1:length(data_list)) {
	id <- names(data_list)[i]
	geneSets[[i]] <- newGeneSet(geneEntrez = data_list[[i]],
				    geneEvidence = "IEA", # Inferred from Electronic Annotation
				    geneSource = "SynGO",
				    ID = id, # diseaseId
				    name = namen[i], # Shortened disease name
				    description = "Synaptic gene ontology",
				    source = myurl,
				    organism = "mouse",
				    internalClassification = c("PL","syngo"),
				    groups = "PL",
				    lastModified = Sys.Date())
}

# Define PLgroup.
PLgroup = newGroup(name = "PL", 
		   description = "Currated synaptic gene ontology from SynGO database.",
		   source = "syngoportal.org")

# Combine as gene collection.
SynGOcollection <- newCollection(dataSets = geneSets, groups = list(PLgroup))

# Save.
myfile <- file.path(rdatdir,paste0(org,"_SynGO_geneSet.RData"))
saveRDS(SynGOcollection,"mouse_SynGO_geneSet.RData")
