#!/usr/bin/env Rscript

# Imports.
suppressPackageStartupMessages({
	library(anRichment)
	library(data.table)
	library(dplyr)
	library(getPPIs)
	library(readxl)
})

# Directories.
here <- getwd()
root <- dirname(here)
rdatdir <- file.path(root,"rdata")
tabsdir <- file.path(root,"tables")
downdir <- file.path(root,"downloads")

# Urls to supplemental data.
urls <- c(S1 = "https://science.sciencemag.org/highwire/filestream/726807/field_highwire_adjunct_files/0/aav8130_Data-S1.xlsx",
	  S2 = "https://science.sciencemag.org/highwire/filestream/726807/field_highwire_adjunct_files/1/aav8130_Data-S2.xlsx",
	  S3 = "https://science.sciencemag.org/highwire/filestream/726807/field_highwire_adjunct_files/2/aav8130_Data-S3.xls",
	  S4 = "https://science.sciencemag.org/highwire/filestream/726807/field_highwire_adjunct_files/3/aav8130_Data-S4.xls",
	  S5 = "https://science.sciencemag.org/highwire/filestream/726807/field_highwire_adjunct_files/4/aav8130_Data-S5.xlsx")

# Description of the data.
descriptions <- list(
		     S1 = "Sample and clinical information for ASD and epilepsy individuals.",
		     S2 = "List of captured nuclei and associated metadata for ASD and epilepsy cohorts.",
		     S3 = "List of cluster-specific and regional gene markers.",
		     S4 = c("List of cell type-specific genes differentially expressed in ASD and epilepsy,",
		            " as well as region-specific and individual-specific gene expression changes in ASD."),
		     S5 = c("Results of whole-exome sequencing analysis of ASD patients.",
		            " The first tab includes high-confidence variants, the second tab",
	    	            " includes variants that are associated with downregulation of",
	                    " corresponding genes in the same ASD patient when compared to control samples;",
			    " the other tabs include unfiltered lists of variants for each individual with less",
			    " confidence in association with ASD, epilepsy or psychiatric disease."))

# Download the data.
sapply(urls,function(x) download.file(x,basename(x)))
myfiles <- basename(urls)

#---------------------------------------------------------------------
## Load the data.
#---------------------------------------------------------------------

## Sheet 1 of S1 has unique formatting. Load it first.
myfile <- file.path(downdir,myfiles[1])
sheets <- excel_sheets(myfile)
S1 <- list()
colNames <- readLines(file.path(downdir,"colnames.txt"))
S1[[1]] <- read_excel(myfile,sheet=1,skip=2,col_names=colNames)
S1 <- c(S1,lapply(sheets[-1],function(x) read_excel(myfile,sheet=x)))
names(S1) <- sheets

# Load the rest of the data. Combine with S1 into a list.
data <- list()
data[[1]] <- S1
data <- c(data, lapply(myfiles[-1],read_sheets))
names(data) <- names(urls)


# Table S4 contains DE genes that are:
# Cell-type specific.
# Region specific (ACC and PFC).
# WGCNA modules.
# Cell-type specific from specific individuals.

# Lets just collect all the genes and see what they are.
S4 <- data$S4
x = S4[[1]]


# Gene identifiers are:
ids <- c("[Gg]ene ID|ENSEMBL ID")
idcols <- sapply(S4,function(x) grep(ids,colnames(x),value=TRUE))

mapply(function(df) df[[col]], S4, idcols)

read_sheets <- function(excel_file) {
	# Function to read all sheets of an excel file.
	# Returns named list of excel sheets.
	sheets <- excel_sheets(excel_file)
	excel_data <- lapply(sheets, function(sheet) read_excel(excel_file,sheet))
	names(excel_data) <- excel_sheets(excel_file)
	return(excel_data)
}
