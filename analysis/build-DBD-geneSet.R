#!/usr/bin/env Rscript

## Building a developmental brain disorder-associated gene set 
# collection for analysis with the AnRichment package.
# Data download from: 
# http://dbd.geisingeradmi.org/#additional-information

# FIXME: script was written to work with "All" dataset, won't work with others. 

## User parameters:
dataset <- "All" # One of c("LOF","Missense","All")

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

# Downoad the raw data.
message("Downloading developmental brain disorder (DBD)-associated gene database... \n")
base_url <- "http://dbd.geisingeradmi.org/downloads"
datasets <- c(LOF = "Full-Gene-Table-Data",
	      Missense = "Full-Missense-Table-Data",
	      All = "DBD-Genes-Full-Data")
myurl <- file.path(base_url,paste0(datasets[dataset],".csv"))
myfile <- paste0(datasets[dataset],".csv")
download.file(myurl,destfile=basename(myurl),quiet=FALSE)
raw_data <- fread(myfile)
unlink(myfile)

# Map human gene symbols to Entrez.
genes <- unique(raw_data$Gene)
entrez <- mapIDs(genes,from="symbol",to="entrez",species="human")
names(entrez) <- genes

# Map missing Entrez IDs by hand.
message("Manually mapping missing ids...")
not_mapped <- genes[is.na(entrez)]
mapped_by_manual_search <- c("HIST1H1E"=3006,"HIST1H4B"=8366,"KIF1BP"=26128)
entrez[not_mapped] <- mapped_by_manual_search[names(entrez[not_mapped])]

# Check.
check <- sum(is.na(entrez)) == 0
if (check)  { message("Successfully mapped all human gene symbols to Entrez!") }
nHsGenes <- length(unique(entrez))

# Add human entrez IDs to data.
data <- tibble::add_column(raw_data,"hsEntrez" = entrez[raw_data$Gene],.after="Gene")

# Map human genes in to their mouse homologs.
message("Mapping human genes to their mouse homologs...\n")
hsEntrez <- data$hsEntrez
msEntrez <- getHomologs(hsEntrez,taxid=10090)
data <- tibble::add_column(data,msEntrez=msEntrez,.after="hsEntrez")
# Remove rows with unmapped genes.
n_out <- sum(is.na(msEntrez))
percent_removed <- round(100*(n_out/length(msEntrez)),2)
message(paste0("Genes without mouse homology: ",n_out,
	      " (",percent_removed," %)."))
data <- data[msEntrez != "NA"] 

# Add concise disorder association column.
disorders <- c("ID/DD","Autism","Epilepsy","ADHD","Schizophrenia","Bipolar Disorder")
idy <- sapply(disorders,function(x) grep(x,colnames(data)))
da <- apply(data[,idy,with=FALSE],1,function(x) {
		    da <- paste(disorders[grep("X",x)],collapse=";")
		    return(da) })
data <- tibble::add_column(data,"disorder_association"=da,.after="msEntrez")

# Separate rows with multiple disorder associations.
data <- tidyr::separate_rows(data,disorder_association,sep=";")

# Fix missing disorder annotation for DYRK1A (Ruaud et al., 2015; ID/DD).
idx <- data$Gene == "DYRK1A" & data$disorder_association == ""
data$disorder_association[idx] <- "ID/DD"

# Drop unnecessary columns.
data[, (disorders):=NULL] 

# Status report.
nMsGenes <- length(unique(data$msEntrez))
nDisorders <- length(unique(data$disorder_association))
message(paste0("Compiled ",nMsGenes," mouse genes associated with ",
	       nDisorders," DBDs!"))

# Write data to file.
myfile <- file.path(tabsdir,paste0("mouse_",datasets[dataset],".csv"))
data.table::fwrite(data,myfile)

# Split into disorder groups.
disorders <- unique(data$disorder_association)
data_list <- data %>% group_by(disorder_association) %>% group_split()
names(data_list) <- disorders

# Check disorder group sizes.
sizes <- sapply(data_list,function(x) length(unique(x$msEntrez)))

# Build gene sets:
geneSets <- list()
for (i in seq_along(data_list)) {
	subdat <- data_list[[i]]
	id <- paste0("DBD-",names(data_list)[i])
	geneSets[[i]] <- newGeneSet(geneEntrez = unique(subdat$msEntrez),
				    geneEvidence = "IEA",
				    geneSource = datasets[dataset],
				    ID = id,
				    name = names(data_list)[i], # disorder name
				    description = "DBD-associated genes.",
				    source = "http://dbd.geisingeradmi.org/",
				    organism = "mouse",
				    internalClassification = c("PL","DBD"),
				    groups = "PL",
				    lastModified = Sys.Date())
} # Ends loop.

# Define gene collection groups.
PLgroup <- newGroup(name = "PL", 
		   description = "DBD-associated genes.",
		   source = "http://dbd.geisingeradmi.org/")

# Combine as gene collection.
DBDcollection <- newCollection(dataSets = geneSets, groups = list(PLgroup))

# Save as RData.
myfile <- file.path(rdatdir,paste0("mouse_",datasets[dataset],"_geneSet.RData"))
saveRDS(DBDcollection,myfile)
