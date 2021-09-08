#!/usr/bin/env Rscript

## using geneLists to perform GSEA

library(geneLists)

# to see all available gene lists:
# geneLists()

library(dplyr)
library(data.table)

# load the ePSD BioID data
data(epsd, package="Uezu2016") # DLG4-BioID sig prots
data(epsd_bioid, package="Uezu2016") # BioID data

# load the SFARI ASD-associated genes
data("sfariGene", package = "geneLists")
data("sfariAnimal", package = "geneLists")

# a geneList is just a named list of Entrez IDs
str(sfariGene)

# combine SFARI gene lists
sfari <- unique(c(sfariGene[["ASD"]], sfariAnimal[["ASD"]]))

# use all genes identified in proteomics experiment as background
all_entrez <- unique(epsd_bioid$Entrez)

# test for enrichment of SFARI genes in DLG4 (epsd) proteome
geneLists::hyperTest(sfari, epsd, background = all_entrez)

# convert entrez IDs to mouse genes using getIDs
genes <- geneLists::getIDs(sfari, "entrez", "symbol", "mouse")
head(genes[which(sfari %in% epsd)])
