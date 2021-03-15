#!/usr/bin/env Rscript

# title: geneLists
# author: twab
# description: creating a geneList

## ---- input 

# the project's root dir (os specific)
root <- "~/projects/soderling-lab/geneLists"

# the name of this script
script <- "064_E3-Ligases" # no file-extension

# the name of the final gene_list object that will be loaded with the data command
short_name <- "e3ligase" 

# the data are from:
data_source <- "https://hpcwebapps.cit.nih.gov/ESBL/Database/E3-ligases/"
data_url <- "https://hpcwebapps.cit.nih.gov/ESBL/Database/E3-ligases/Definite%20Ligase%20List.xlsx"

# this source includes data from Li2008:
doi <- "https://doi.org/10.1371/journal.pone.0001487"

# FIXME: geneLists should be associated with reference and possibly abstract
# info! Using doi or pubmed there should be automatic ways to get abstract...
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


# for some misc functions in this package
# * documentDataset
# * write_gmt
devtools::load_all(root, quiet=TRUE)


## ---- download the data

# the downloads dir is not tracked by git
stopifnot(dir.exists(file.path(root,"downloads")))

# download the data
myfile <- file.path(root,"downloads","ligases.xlsx")
download.file(data_url, destfile=myfile)

# read into R
df <- read_excel(myfile)

# we need the human gene map, install it if required
#BiocManager::install("org.Hs.eg.db")

# apparently these arent refseq...
#entrez <- getIDs(genes,from="refseq",to="entrez","human")

# and Swiss-Prot is actually uniprot, yes gene identifiers are confusing!

# we can map the uniprot IDs to entrez
uniprot <- df$"Swiss-Prot"

# map to uniprot
df$hsEntrez <- getIDs(uniprot,"uniprot","entrez","human")

# how many were not mapped
sum(is.na(df$hsEntrez))

# drop 8 NA
subdf <- df %>% filter(!is.na(hsEntrez))

# map to mouse with getHomologs
subdf <- subdf %>% mutate(msEntrez = geneLists::getHomologs(hsEntrez,"mouse"))

# 32 not mapped to homologs
sum(is.na(subdf$msEntrez))

head(subdf)

# create gene_list excluding NA
is_na <- is.na(subdf$msEntrez)
groups <- subdf$Domain[!is_na]
gene_list <- split(subdf$msEntrez[!is_na],groups)

# check...
subdf %>% filter(`Gene Symbol` == "UBE3A")

# categories of ligases
l <- sapply(gene_list, length) 
data.table(group=names(l), k=l) %>% knitr::kable()

x = unlist(gene_list)
length(x)


## ---- save

# write gmt file, saved in root/gmt
#write_gmt(gene_list, data_source, file.path(root,"gmt",script))

# generate simple dataset documentation in root/R 
#   and generate rda object in root/data
#documentDataset(file.path(root,"gmt",script), short_name, 
#		file.path(root,"R"), file.path(root,"data"))

# now push changes to git and you can use
#data(e3ligase)
