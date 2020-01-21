#!/usr/bin/env Rscript

# Interneuron iBioID enrichment analysis. 

library(readxl)
library(getPPIs)

# Directories.
here <- getwd()
root <- dirname(here)
funcdir <- file.path(root,"R")
rdatdir <- file.path(root,"rdata")
downdir <- file.path(root,"downloads")

# Source functions.
invisible(sapply(list.files(funcdir,full.names=TRUE),source))

# Load top iBioID proteins.
myfile <- file.path(downdir,"topSig.xlsx")
data <- list()
for (sheet in excel_sheets(myfile)){
	data[[sheet]] <- read_excel(myfile,sheet)
}

# Insure we only have significantly enriched genes.
cutoff <- 0
alpha <- 0.05
data <- lapply(data,function(x) filter(x,FDR < alpha & logFC > cutoff))

# Number of significant genes.
ngenes <- sapply(data,dim)[1,]

# Collect genes.
uniprot <- lapply(data,function(x) x$Accession)

# Map uniprot to entrez.
entrez <- lapply(uniprot,function(x) {
			 mapIDs(x,from="uniprot",to="entrez",species="mouse")
})

# Load custom gene lists.
myfiles <- list.files(rdatdir,full.names=TRUE)
geneSets <- lapply(myfiles,readRDS)
names(geneSets) <- gsub(".RData","",basename(myfiles))

# Add Mouse GO database.
musGO <- buildGOcollection(organism="mouse")
geneSets[["mouse_GO"]] <- musGO

# Perform gene enrichment analysis for all gene sets.
for (i in seq_along(geneSets)){
	if (i == 1) { pbar <- txtProgressBar(min=1,max=length(geneSets),style=3) } 
	setTxtProgressBar(pbar,i)
	results[[i]] <- gene_list_enrichment(entrez, geneSets[[i]])
	if (i == length(geneSets)) { close(pbar); message("\n") }
}
names(results) <- names(geneSets)
all_results <- unlist(results,recursive=FALSE)

# Shorten names... max length is 31 characters.
new_names <- gsub("_geneSet","",gsub("mouse_","",names(all_results)))
new_names <- gsub("_All_Disease_Genes","",new_names)
names(all_results) <- new_names

# Write to file.
write_excel(all_results,"Interneuron_GeneSet_Enrichment.xlsx")
