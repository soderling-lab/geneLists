#!/usr/bin/env Rscript

## User parameters:

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
myfiles <- list.files(rdatdir,pattern="geneSet",full.names=TRUE)
geneSets <- lapply(myfiles,readRDS)
names(geneSets) <- tools::file_path_sans_ext(basename(myfiles))

# Combine.
combinedCollection <- do.call(mergeCollections,geneSets)
