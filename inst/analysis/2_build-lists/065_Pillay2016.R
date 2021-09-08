#!/usr/bin/env Rscript

# title: geneLists
# author: twab
# description: compiling genes involved in AAV transduction by Pillay2016

## ---- input 

# the project's root dir (os specific)
root <- "~/projects/soderling-lab/geneLists"
devtools::load_all(root)

# the name of this script
script <- "065_Pillay2016" # no file-extension
short_name <- "Pillay2016" # name for rda object

# the data are from:
doi <- "10.1038/nature16465"
data_source <- "https://pubmed.ncbi.nlm.nih.gov/26814968/"

# todo:
# keywords!
# class with methods:
# * search for a genelist by keywords: geneLists(keywords="AAV")
# * cite a genelist: Pillay2016.cite or ref(Pillay2016)
# * keep track of source organism details: Pillay2016.organism
# * always keep source organism genes! 
# * interconvert to another organism: Pillay2016.convert(to="mouse")


## ---- prepare the renv

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(data.table)
})


## ---- work with reference

# this function should be internal
bib <- doi2bibentry(doi)

print(bib,style="bibtex")


## ---- load the data

#download.file(data_url, destfile=myfile)
myfile <- list.files(file.path(root,"downloads"),pattern="Pillay2016-SI*", full.name=TRUE)
df <- data.table::fread(myfile)

# fix bad gene symbols 
myfile <- list.files(file.path(root,"downloads"),pattern="Pillay2016-map*", full.name=TRUE)
gene_map <- data.table::fread(myfile)
idx <- match(gene_map$symbol, df$Symbol)
df$Symbol[idx] <- gene_map$gene

# map to entrez
symbols  <- df$Symbol
entrez <- getIDs(symbols,from="symbol",to="entrez","human")
df$entrez <- entrez

# check mapping
# NOTE: there is one gene (lncRNA) that is not mapped to entrez
df %>% filter(is.na(entrez)) %>% dplyr::select(Symbol,entrez)

# drop NA
subdf <- df %>% filter(!is.na(entrez))

# map to mouse
#subdf <- subdf %>% mutate(msEntrez = geneLists::getHomologs(entrez,"human","mouse"))

# n not mapped
#sum(is.na(subdf$msEntrez))

# the authors report 46 significant hits (alpha = 0.001)
idx <- subdf$`q-value` < 0.001

gene_list <- geneList("Pillay2016", subdf$entrez[idx], "human", bib)

# an abstract may not always be include...
#pub_data <- rcrossref::cr_works(dois = doi)[["data"]]


## ---- save

# write gmt file, saved in root/gmt
write_gmt(gene_list, data_source, file.path(root,"gmt",script))

# generate simple dataset documentation in root/R 
#   and generate rda object in root/data
# FIXME: this is were the problem arrises. The gene list is written to file as
# GMT and then read and saved. The document dataset function adds dataset name,
# url to template file...
documentDataset(file.path(root,"gmt",script), short_name, 
		file.path(root,"R"), file.path(root,"data"))
