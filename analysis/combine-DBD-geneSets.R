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
patterns <- c("SFARI-Gene_genes","SFARI-Gene_animal","DDG2P",
	      "DisGeneNet_All_Disease_Genes","URMC","DBD-Genes")
myfiles <- sapply(patterns,function(x) list.files(rdatdir,x,full.names=TRUE))
geneSets <- lapply(myfiles,readRDS)
names(geneSets) <- tools::file_path_sans_ext(basename(myfiles))

# Combine.
combinedCollection <- do.call(mergeCollections,geneSets)

# Save.
myfile <- file.path(rdatdir,paste0("mouse_Combined_DBD_geneSets.RData"))
saveRDS(combinedCollection,myfile)
