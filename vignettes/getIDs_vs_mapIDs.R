#!/usr/bin/env Rscript

# title: geneLists
# author: twab
# description: compare usage of getIDs versus mapIDs

library(dplyr)
library(data.table)

library(geneLists) # for getIDs and mapIDs

library(SwipProteomics)

data(partition)
# uniprot IDS
head(partition)

# create data.table
df = data.table(prot = names(partition))

# use getIDs to map Uniprot to gene symbols. 
# note: from and to arguments are not case sensitive
# note: but, watch out for symbol versus symbols plural!
df = df %>% mutate(symbol = getIDs(prot, from="uniprot",to="symbol", "mouse"))

# getIDs using org.mm.eg.db for mouse
head(df)

# UniprotIDs are dynamic and hard to track!
sum(is.na(df$symbol))

# use mapIDs with a user defined gene_map

data(swip_gene_map) # == gene_map
# ^note: its probably a better habit to save data with the same name of the R object!

head(gene_map) 
# ^note: you might have to create this object by hand using multiple methods to successfuly
# map all uniprot to another identifier

df = df %>% mutate(symbol = mapIDs(prot,"uniprot","symbol",gene_map))

# no missing
sum(is.na(df$symbol))
