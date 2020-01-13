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
myfiles <- list.files(rdatdir,pattern=c(org,"geneSet"),full.names=TRUE)
geneSets <- lapply(myfiles,readRDS)
names(geneSets) <- tools::file_path_sans_ext(basename(myfiles))

# Combine.
combinedCollection <- do.call(mergeCollections,geneSets)

# Save.
myfile <- file.path(rdatdir,paste0(org,"_Combined_geneSets.RData"))
saveRDS(combinedCollection,myfile)
