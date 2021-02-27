#!/usr/bin/env Rscript

# title: geneLists
# author: twab
# description: using get homologs

library(geneLists)

# ?getHomologs

install.packages("org.Ce.eg.db")

# c elegans lin23 uniprot
lin23 = getIDs("Q09990","uniprot","entrez","worm")

# get mouse homolog
getHomologs(lin23,"mouse")

# get human homolog
getHomologs(lin23,"human")

# map mouse homolog to gene symbol
getIDs("103583","entrez","symbol","mouse")
