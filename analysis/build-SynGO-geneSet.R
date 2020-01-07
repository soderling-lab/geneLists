#!/usr/bin/env Rscript

## Building Synaptic GO collection for GO analysis with the 
# AnRichment package.
#  SynGO downloaded from: https://www.syngoportal.org/ 

suppressPackageStartupMessages({
	library(anRichment)
})

# Load SynGO.
myfile <- "SynGO_Pathways.RData"
data_list <- do.call(c,readRDS(myfile))

# Arbitrary names.
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
				    source = "https://syngoportal.org/data/download.php?file=SynGO_bulk_download_release_20180731.zip",
				    organism = "mouse",
				    internalClassification = c("PL","syngo"),
				    groups = "PL",
				    lastModified = "2020-01-03")
}

# Define PLgroup.
PLgroup = newGroup(name = "PL", 
		   description = "Currated synaptic gene ontology from SynGO database.",
		   source = "syngoportal.org")
# Combine as gene collection.
SynGOcollection <- newCollection(dataSets = geneSets, groups = list(PLgroup))

# Save.
saveRDS(SynGOcollection,"SynGOcollection.RDS")
