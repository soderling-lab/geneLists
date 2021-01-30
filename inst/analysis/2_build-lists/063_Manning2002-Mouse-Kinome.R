#!/usr/bin/env Rscript

## ---- input 

root <- "~/projects/geneLists"
script <- "063_Manning2002-Human-Kinome" # no file-extension
short_name <- "msKinome"
data_source <- "https://pubmed.ncbi.nlm.nih.gov/12471243/"
data_url <- "http://kinase.com/static/colt/data/human/kinome/tables/Kincat_Hsap.08.02.xls"

## The best article identifier is its DOI!
# * DOI: 10.1126/science.1075762

## ---- documentation


## ---- prepare the renv


# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

devtools::load_all(root, quiet=TRUE)

data(hsKinome)


## ---- load the data

df <- data.table(entrez = unlist(hsKinome,use.names=F),
		 group = unlist(sapply(names(hsKinome),function(x) rep(x,length(hsKinome[[x]])))))

hs_entrez <- df$entrez


## ---- map to mouse

# get mouse homologs
mus_entrez <- geneLists::getHomologs(hs_entrez,species="mouse")

# create gene_list
is_na <- is.na(mus_entrez)
groups <- df$group[!is_na]
gene_list <- split(mus_entrez[!is_na],groups)

# drop na, get unique genes
gene_list <- lapply(gene_list, function(x) unique(x[!is.na(x)]))

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
