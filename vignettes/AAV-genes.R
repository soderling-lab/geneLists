#!/usr/bin/env Rscript

# title:
# author:
# description:

## using geneLists to perform GSEA

devtools::load_all()
#data("Pillay2016", package = "geneLists")

data("swip_partition", package = "SwipProteomics")
data("swip_gene_map", package = "SwipProteomics")

entrez <- geneLists::mapIDs(names(partition),from="uniprot",to="entrez",gene_map)
modules <- split(entrez,partition)
names(modules) <- paste0("M",names(modules))

goi <- geneLists::getHomologs(unlist(Pillay2016),species="mouse")
goi <- goi[!is.na(goi)]

# test for enrichment of AAV transduction genes in swip spatial proteome
geneLists::hyperTest(goi, modules[["M22"]], background = c(goi,entrez))
