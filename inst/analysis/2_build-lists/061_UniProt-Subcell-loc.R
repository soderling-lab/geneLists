#!/usr/bin/env Rscript

## ---- input 
root <- "~/projects/geneLists"
script <- "061_UniProt-Subcell-loc"
short_name <- "uniprotSubcell"
data_source <- "uniprot.org"

## ---- prepare the renv

devtools::load_all(root, quiet=TRUE)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})


# downloaded Jan 20th 2021
myfile <- file.path(root,"downloads","uniprot-subcell-01202021.tab.gz")
stopifnot(file.exists(myfile))

# if you have R.utils, you can read .gz directly!
df <- fread(myfile) %>% 
	# munge
	dplyr::rename(Accession = "Entry") %>%
	dplyr::rename(CC = "Subcellular location [CC]") %>%
	select(Accession,CC) %>%
	tidyr::separate_rows(CC,sep="\\;") %>%
	filter(CC != "") %>%
	mutate(CC = gsub("SUBCELLULAR LOCATION: ","", CC)) %>%
	mutate(CC = trimws(CC))

# need to delimit two cases:
# multiple entries seperated by "."
# do not delimit ". Note="
subdf = df %>% dplyr::filter(grepl("Note=",CC))

all(grepl("Note",subdf$CC))

# the notes look really useful, but we cannot keep this information
subdf <- subdf %>% tidyr::separate(CC,into=c("CC","Note"),sep="Note=")

# now we can seperate CC
subdf <- subdf %>% 
	tidyr::separate_rows(CC,sep="\\. ") %>%
	filter(CC!="")

# there are multiple types of entries (all associated with notes)
# * just loc: 'Cytoplasm'
# * loc and evidence 'Cytopasm {}

# there are multiple types of evidence annot (ECO) 
# https://www.uniprot.org/help/evidences
subdf2 = df %>% 
	dplyr::filter(!grepl("Note=",CC)) %>% 
	mutate(Note="") %>%
	tidyr::separate_rows(CC,sep="\\. ") %>%
	filter(CC!="")

# combine and final clean-up
res <- rbind(subdf,subdf2) %>% 
	mutate(CC = gsub("\\.","",CC)) # rm trailing "."


# okay, now lets group things together
part1 <- res %>% dplyr::filter(grepl("\\{",CC)) # evidence
part2 <- res %>% dplyr::filter(!grepl("\\{",CC)) # no evidence

# we want to parse things into nice bins


# to account for cases like the following row:
# "[Amyloid-beta protein 42]: Cell surface"
# we can get subcell annot from the part after []
subpart1 <- part1 %>% 
	filter(grepl("\\[",CC)) %>%
	tidyr::separate(CC,into=c("Isoform","CC"),sep="\\[.*\\]") %>%
	mutate(CC=gsub("\\:","",CC)) %>%
	mutate(CC=trimws(CC)) %>%
	# we loose isoform specific info... could be appended to note
	select(-Isoform) %>%
	tidyr::separate(CC,into=c("CC","Evidence"),sep="\\{") %>%
	mutate(CC = trimws(CC)) %>%
	mutate(Evidence=gsub("\\}","",Note))

subpart2 <- part1 %>%
	filter(!grepl("\\[",CC)) %>%
	tidyr::separate(CC,into=c("CC","Evidence"),sep="\\{") %>%
	mutate(CC = trimws(CC)) %>%
	filter(CC != "") %>%
	mutate(Evidence=gsub("\\}","",Note))

# check how many locs
unique(subpart1$CC)
unique(subpart2$CC)

# now for the other half
# part2

# NOTE: evidence can be buried in note!
# we ignore this additional complexity for now

subpart3 <- part2 %>% 
	filter(grepl("\\[",CC)) %>%
	tidyr::separate(CC,into=c("Isoform","CC"),sep="\\[.*\\]") %>%
	mutate(CC=gsub("\\:","",CC)) %>%
	mutate(CC=trimws(CC)) %>%
	# we loose isoform specific info... could be appended to note
	mutate(Evidence=NA) %>%
	select(-Isoform) %>%
	select(Accession,CC,Evidence,Note)

subpart4 <- part2 %>% 
	filter(!grepl("\\[",CC)) %>%
	mutate(Evidence=NA) %>%
	select(Accession,CC,Evidence,Note)


# whew!
uniprot_subcell <- rbind(subpart1,subpart2,subpart3,subpart4)

# note there are some strange CC anno coming from subpart4...
# these result in groups of size == 1


uniprot <- uniprot_subcell$Accession

# n missing is about the same for both approaches

#entrez <- geneLists::getIDs(uniprot,"uniprot","entrez","mouse")
#sum(is.na(entrez))

entrez <- geneLists::mapIDs(uniprot,"Accession","Entrez",uniprot_map)
#sum(entrez=="")

uniprot_subcell <- uniprot_subcell %>% 
	mutate(Entrez = entrez) %>% 
	filter(Entrez != "") %>%
	filter(!is.na(Entrez))


## ---- create gene list

gene_list <- split(uniprot_subcell$Entrez,uniprot_subcell$CC)

too_small <- names(which(sapply(gene_list,length) < 3))

gene_list <- gene_list[!(names(gene_list) %in% too_small)]

message("Compiled n subcellular locations (pathways): ", length(gene_list))

df <- data.table(Pathway=names(gene_list),n=sapply(gene_list,length))
df %>% head() %>% knitr::kable()

#namen <- sample(names(gene_list),1)
#path <- gene_list[[namen]]
#data.table(Pathway=namen,Genes=length(path))
#head(getIDs(path,"entrez","symbol","mouse"))

myfile <- file.path(root,"gmt",script)
write_gmt(gene_list,data_source,myfile)
documentDataset(myfile,short_name,file.path(root,"R"),file.path(root,"data"))
