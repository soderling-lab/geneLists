#!/usr/bin/env Rscript

## Annotating CORUM protein complexes with predicted subcellular localization.

# Load renv.
renv::load(getrd())

# Load geneLists library and other dependencies.
library(geneLists)
library(dplyr) # For manipulating the data.
library(data.table) # For working with tables.

# Load corum dataset.
data(corum)
head(corum)

# Load subcellular predictions:
data(orre2019concensus) # Concensus predictions from Orre et al., 2019.
data(itzhak2016predictions) # Subcellular predictions from Itzhak et al., 2016.
data(itzhak2017predictions) # Subcellular predictions from Itzhak et al., 2017.
data(lopitDCpred) # Subcellular predictions from Geledaki et al., 2019.

# Let's combine these datasets.
markers <- list(orre19 = orre2019concensus,
		itzhak16 = itzhak2016predictions,
		itzhak17 = itzhak2017predictions,
		geledaki19 = lopitDCpred)

# Combine the redundant categories.
# Capitalize names of each gene list.
capitalize_names <- function(x) { 
	names(x) <- toupper(names(x)) 
	return(x)
}
markers <- lapply(markers,capitalize_names)

# Convert each list of genes to a data.frame.
geneList2df <- function(gene_list,class.name="class"){
	entrez <- unlist(gene_list,use.names=FALSE)
	names(entrez) <- rep(names(gene_list), times = sapply(gene_list,length))
	df <- data.frame(class=names(entrez),entrez)
	colnames(df)[1] <- class.name
	return(df)
}
df_list <- lapply(markers,geneList2df)

# Combine these dataframes into a single df.
df <- bind_rows(df_list,.id="column_label")
colnames(df)[which(colnames(df)=="column_label")] <- "study"

# Clean-up class names.
df$class <- gsub("ER_HIGH_CURVATURE|ER$","ENDOPLASMIC RETICULUM",df$class)
df$class <- gsub("UNCLASSIFIED","UNKNOWN",df$class)
df$class <- gsub("GA$","GOLGI",df$class)
df$class <- gsub("PM$","PLASMA MEMBRANE",df$class)
df$class <- gsub("NUCLEAR PORE COMPLEX\\/NUCLEAR$","NUCLEUS",df$class)
df$class <- gsub("NUCLEAR PORE COMPLEX$","NUCLEUS",df$class)
df$class <- gsub("NUCLEUS\\/CHROMATIN$","NUCLEUS",df$class)
df$class <- gsub("NUCLEAR$","NUCLEUS",df$class)
df$class <- gsub("MITOCHONDRION$","MITOCHONDRIA",df$class)
df$class <- gsub("ERGIC\\/CISGOLGI$","GOLGI",df$class)
df$class <- gsub("ACTIN BINDING PROTIENS$","ACTIN BINDING",df$class)
marker_df <- df
all_classes <- unique(marker_df$class)
print(all_classes)

# Summarize number of proteins per class:
marker_df %>% group_by(class) %>% 
	summarize(n = length(unique(entrez))) %>%
	arrange(desc(n))

# Coerce corum list to a data.frame.
corum_df <- geneList2df(corum,class.name="complex")
head(corum_df)

# Combine corum data with subcellular predictions.
data <- full_join(corum_df,marker_df, by="entrez") %>% as.data.table()
data <- data %>% filter(!is.na(complex))
data_list <- data %>% group_by(complex) %>% group_split()
names(data_list) <- sapply(data_list,function(x) unique(x$complex))

# Predict a complexes localization based on the most frequent subcellular
# annotation.
get_max_pred <- function(df) {
	prediction <- names(which(table(df$class) == max(table(df$class))))
	return(paste(prediction,collapse="|"))
}

# For every complex, predict its subcellular localization.
pred <- vector(mode="character",length(data_list))
names(pred) <- names(data_list)
for (i in 1:length(data_list)){
	df <- data_list[[i]]
	# Use try catch to handle errors raised when none of the complex
	# proteins have a predicted subcellular localization.
	pred[i] <- tryCatch(get_max_pred(df), warning=function(w) { return("UNCLASSIFIED") })
}
pred_loc <- data.frame(complex=names(pred),prediction=pred)

# Combine with other corum data.
df_a <- corum_df %>% group_by(complex) %>%
	summarize(n=length(unique(entrez)),
		  entrez=paste(unique(entrez),collapse="; "))
results <- left_join(pred_loc,df_a,by="complex")

# Write to file.
fwrite(results,"complex_predicted_localization.csv")
