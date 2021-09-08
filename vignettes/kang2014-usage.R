#!/usr/bin/env Rscript

# analyze kang2014 proteome

## ---- imports 

library(dplyr)
library(data.table)

library(geneLists)


## ---- load the data

data(ciPSD, package="geneLists")
data(uniprotSubcell, package="geneLists")


## ---- main

kang2014 = ciPSD[["Kang2014"]]

# goi
kang2014

# background
back <- unique(c(kang2014, unlist(uniprotSubcell)))

# screen for uniprot subcell enrichment
results <- list()
for (i in c(1:length(uniprotSubcell))) {
	path <- uniprotSubcell[[i]]
	namen <- names(uniprotSubcell)[i]
	res <- as.data.frame(t(hyperTest(path, kang2014, background = back))) %>% mutate(Pathway=namen)
	results[[i]] <- res
}

#df = do.call(rbind,results) %>% arrange(desc(`Fold enrichment`))
#fwrite(df, "foo.csv")

