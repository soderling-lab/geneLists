#!/usr/bin/env Rscript

# Imports.
suppressPackageStartupMessages({
	library(anRichment)
})

# Directories.
here <- getwd()
root <- dirname(here)
rdatdir <- file.path(root,"rdata")
tabsdir <- file.path(root,"tables")
downdir <- file.path(root,"downloads")

# Load geneSets.
datasets <- c("SFARI-Gene_genes",
	      "SFARI-Gene_animal",
	      "DDG2P",
	      "DisGeneNet_All_Disease_Genes",
	      "URMC",
	      "DBD-Genes",
	      "Satterstrom",
	      "Sanders_ASD_geneSet")
myfiles <- sapply(datasets,function(x) list.files(rdatdir,x,full.names=TRUE))
geneSets <- lapply(myfiles,readRDS)
names(geneSets) <- tools::file_path_sans_ext(basename(myfiles))

# Some summary stats.
data_list <- lapply(geneSets,function(x) x$dataSets[[1]]$data)
df <- do.call(rbind,data_list)
n_associations <- dim(df)[1]
n_genes <- length(unique(df$Entrez))
n_datasets <- length(geneSets)
message(paste("Compiled",n_associations,"DBD-gene associations from",
	      n_datasets, "databases cooresponding to",
	      n_genes, "unique DBD-associated mouse genes."))

# Combine.
combinedCollection <- do.call(mergeCollections,geneSets)

# Save.
myfile <- file.path(rdatdir,paste0("mouse_Combined_DBD_geneSets.RData"))
saveRDS(combinedCollection,myfile)
