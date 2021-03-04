#!/usr/bin/env Rscript

# title: geneLists
# author: twab
# description: creating a geneList

## ---- input 

# the project's root dir (os specific)
root <- "~/projects/soderling-lab/geneLists"

# the name of this script
script <- "064_E3-Ligases" # no file-extension
short_name <- "e3ligase"

# the data are from:
data_source <- "https://hpcwebapps.cit.nih.gov/ESBL/Database/E3-ligases/"
data_url <- "https://hpcwebapps.cit.nih.gov/ESBL/Database/E3-ligases/Definite%20Ligase%20List.xlsx"

# this includes data from Li2008:
doi <- "https://doi.org/10.1371/journal.pone.0001487"

#@article{Li_2008,
#	doi = {10.1371/journal.pone.0001487},
#	url = {https://doi.org/10.1371%2Fjournal.pone.0001487},
#	year = 2008,
#	month = {jan},
#	publisher = {Public Library of Science ({PLoS})},
#	volume = {3},
#	number = {1},
#	pages = {e1487},
#	author = {Wei Li and Mario H. Bengtson and Axel Ulbrich and Akio Matsuda and Venkateshwar A. Reddy and Anthony Orth and Sumit K. Chanda and Serge Batalov and Claudio A. P. Joazeiro},
#	editor = {Hidde Ploegh},
#	title = {Genome-Wide and Functional Annotation of Human E3 Ubiquitin Ligases Identifies {MULAN}, a Mitochondrial E3 that Regulates the Organelle{\textquotesingle}s Dynamics and Signaling},
#	journal = {{PLoS} {ONE}}
#}

## ---- prepare the renv

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(data.table)
})

# for gene id mapping
library(geneLists)

# for some misc functions in this package
devtools::load_all(root, quiet=TRUE)


## ---- load the data

stopifnot(dir.exists(file.path(root,"downloads")))

myfile <- file.path(root,"downloads","ligases.xlsx")

download.file(data_url, destfile=myfile)

df <- read_excel(myfile)

# these are human refseq identifiers, lets map to mouse
# first map to human entrez, then map to mouse homologs
genes <- df$RefSeq

# we need the human gene map
#BiocManager::install("org.Hs.eg.db")
entrez <- getIDs(genes,from="refseq",to="entrez","human")

# apparently these arent refseq...

# we can map the uniprot IDs instead or gene symbols...
uniprot <- df$"Swiss-Prot"

df$entrez <- getIDs(uniprot,"uniprot","entrez","human")

# drop 8 NA
subdf <- df %>% filter(!is.na(entrez))

# map to mouse
subdf <- subdf %>% mutate(msEntrez = geneLists::getHomologs(entrez,"human","mouse"))

sum(is.na(subdf$msEntrez))
# 2 not mapped

# create gene_list
is_na <- is.na(subdf$msEntrez)
groups <- subdf$Domain[!is_na]

gene_list <- split(subdf$msEntrez[!is_na],groups)

df %>% filter(`Gene Symbol` == "UBE3A")

# categories of ligases
names(gene_list)

# summarize
l <- sapply(gene_list,length) 
data.table(group=names(l),k=l) %>% knitr::kable()


## ---- save

# write gmt file, saved in root/gmt
write_gmt(gene_list, data_source, file.path(root,"gmt",script))

# generate simple dataset documentation in root/R 
#   and generate rda object in root/data
documentDataset(file.path(root,"gmt",script), short_name, 
		file.path(root,"R"), file.path(root,"data"))

# now push changes to git and you can use
#data(e3ligase)
