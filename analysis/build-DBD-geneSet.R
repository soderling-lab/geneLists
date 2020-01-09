#!/usr/bin/env Rscript

## Building a developmental brain disease-associated gene set 
# collection for analysis with the AnRichment package.
# Data download from: 
# http://dbd.geisingeradmi.org/#additional-information

# FIXME: script was written to work with "All" dataset, won't work with others. 

## User parameters:
dataset <- "All" # One of c("LOF","Missense","All")
map2mouse <- TRUE
save_data <- TRUE

# Imports.
suppressPackageStartupMessages({
	library(anRichment)
	library(data.table)
	library(dplyr)
	library(getPPIs)
})

# Directories.
here <- getwd()
root <- dirname(here)
rdatdir <- file.path(root,"rdata")
tabsdir <- file.path(root,"tables")
downdir <- file.path(root,"downloads")

# Which organism, human or mouse?
if (map2mouse) { org <- "mouse" } else { org <- "human" }

# Downoad the data.
base_url <- "http://dbd.geisingeradmi.org/downloads"
datasets <- c(LOF = "Full-Gene-Table-Data",
	      Missense = "Full-Missense-Table-Data",
	      All = "DBD-Genes-Full-Data")
myurl <- file.path(base_url,paste0(datasets[dataset],".csv"))
myfile <- paste0(datasets[dataset],".csv")
download.file(myurl,destfile=basename(myurl),quiet=FALSE)
data <- fread(myfile)
unlink(myfile)

# Map human gene symbols to Entrez.
genes <- unique(data$Gene)
entrez <- mapIDs(genes,from="symbol",to="entrez",species="human")
names(entrez) <- genes

# Map missing Entrez IDs by hand.
not_mapped <- genes[is.na(entrez)]
mapped_by_manual_search <- c("HIST1H1E"=3006,"HIST1H4B"=8366,"KIF1BP"=26128)
entrez[not_mapped] <- mapped_by_manual_search[names(entrez[not_mapped])]

# Check.
check <- sum(is.na(entrez)) == 0
if (check)  { message("Successfully mapped all human gene symbols to Entrez!") }

# Add entrez IDs to data.
data <- tibble::add_column(data,"Entrez" = entrez[data$Gene],.after=1)

# Map human genes in to their mouse homologs.
if (map2mouse) {
	message("Mapping human genes to their mouse homologs...\n")
	hsEntrez <- data$Entrez
	msEntrez <- getHomologs(hsEntrez,taxid=10090)
	data <- tibble::add_column(data,osEntrez=msEntrez,.after=2)
	# Remove rows with unmapped genes.
	n_out <- sum(is.na(msEntrez))
	percent_removed <- round(100*(n_out/length(msEntrez)),2)
	message(paste0("Genes without mouse homology: ",n_out,
		      " (",percent_removed," %)."))
	data <- data[osEntrez != "NA"] 
} else {
	hsEntrez <- data$entrez_id
	data <- tibble::add_column(data,osEntrez=osEntrez,.after=2)
}

# Add concise disease association column.
da <- apply(data[,c(1:9)],1,function(x) {
		    da <- paste(colnames(data)[grep("X",x)],collapse=";")
		    return(da) })
data <- tibble::add_column(data,"disease_association"=da,.after=4)

# Fix missing disease annotation for DYRK1A (Ruaud et al., 2015; ID/DD).
idx <- data$Gene == "DYRK1A" & data$disease_association == ""
data$disease_association[idx] <- "ID/DD"

# Drop unnecessary columns.
cols_out <- colnames(data)[c(6:9)]
data[, (cols_out):=NULL] 

# Note some genes are annotated as disease_association:"Gene".
# These genes arise from a number of studies linking them to ID/DD,
# autism, epilepsy, and other congenital defects. We will just call 
# these "unannoated" DBD genes.
data$disease_association[data$disease_association == "Gene"] <- "unannotated"

# Write data to file.
if (save_data) {
	myfile <- file.path(tabsdir,paste0(org,"_",datasets[dataset]))
	data.table::fwrite(data,myfile)
}

# Split into disease groups.
diseases <- unique(data$disease_association)
nDiseases <- length(diseases)
data_list <- data %>% group_by(disease_association) %>% group_split()
names(data_list) <- diseases

# Check disease group sizes.
sizes <- sapply(data_list,function(x) length(unique(x$osEntrez)))

# Build gene sets:
geneSets <- list()
for (i in seq_along(data_list)) {
	subdat <- data_list[[i]]
	id <- paste0("DBD-",names(data_list)[i])
	geneSets[[i]] <- newGeneSet(geneEntrez = unique(subdat$osEntrez),
				    geneEvidence = "IEA",
				    geneSource = datasets[dataset],
				    ID = id,
				    name = names(data_list)[i], # disease name
				    description = "DBD-associated genes.",
				    source = "http://dbd.geisingeradmi.org/",
				    organism = org,
				    internalClassification = c("PL","DBD"),
				    groups = "PL",
				    lastModified = Sys.Date())
}

# Define gene collection groups.
PLgroup <- newGroup(name = "PL", 
		   description = "DBD-associated genes.",
		   source = "http://dbd.geisingeradmi.org/")

# Combine as gene collection.
DBDcollection <- newCollection(dataSets = geneSets, groups = list(PLgroup))

# Save as RData.
myfile <- file.path(rdatdir,paste0(org,"_",datasets[dataset],".RData"))
saveRDS(DBDcollection,myfile)
