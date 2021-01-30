#!/usr/bin/env Rscript

## Perform hypergeometric test for enrichment.

# Imports.
library(geneLists)

# Load some gene lists.
data(wang2017Epilepsy)
data(synaptosome)
data(sharma2015brain)

# Test enrichment of epilepsy genes in iPSD
epilepsy_genes <- unique(unlist(wang2017Epilepsy))
pathway <- iPSD
brain <- sharma2015brain$All_Brain

# hyperTest - perform hypergeometric test.
hyperTest(epilepsy_genes,iPSD,background=brain_proteome)
