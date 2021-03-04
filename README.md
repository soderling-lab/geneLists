# geneLists

_geneLists_ is an __R__ package containing a collection of
gene lists for simple, reproducible Gene Set Enrichment Analysis (GSEA).


## Description
A major focus of this repository is the collection of __synaptic proteome__
genes as well as genes that are implicated in human __brain disorders__.

Gene lists are stored in the Broad Institute's [GMT
format](https://bit.ly/38gtRLW).  These can be downloaded directly from the
`datasets/` directory, or accessed in __R__ with the `data()` command.  For
example, load the [SFARI](https://gene.sfari.org/) autism candidate gene dataset
with `data(sfariGene)`.

Gene lists are are collected from the literature or online databases. Gene
identifiers are mapped to stable, unique Entrez IDs. Often, it is
necessary to map human genes to their homlogous mouse genes. This is done using
the Homologene database and the `getHomologs` function.

## Installation
Insure you have installed `AnnotationDbi` beforing installing `geneLists`.
For example in __R__, download `geneLists` from GitHub using the devtools package:

```R
# Install from github
devtools::install_github("soderling-lab/geneLists")
```

## Usage

```R
library(geneLists)

# See all available datasets
geneLists()

# Load a dataset
data(iPSD) # Uezu2016 iPSD genes

# converting between identifiers
gphn_proteome <- iPSD[["Gphn"]]
uniprot <- getIDs(gphn_proteome, from="entrez", to="uniprot", species="mouse")

# mapping genes using a given gene map
data(uniprot_map)
mapIDs(uniprot,"Accession", "Entrez", uniprot_map)

# NOTE: be careful to not confuse getIDs (uses org.##.eg.db) and
# mapIDs (you must provide a gene_map)

```

## Datasets
For additional details about each dataset, see the [README](./datasets/README.md)
in the `datasets/` directory. The source code used to compile each gene list is
in `inst/analysis`.

```R
# to see all scripts in inst/analaysis/2_build-lists:
list.files(system.file("analysis/2_build-lists", package="geneLists"))
```

The scripts in `inst/analysis` record how each gene list was created. These
can be used as examples to show you can download a dataset and save it as a GMT
formatted file and gene_list R object. See the `tutorials/064_E3-Ligases.R`
script for a recent example on how to create a gene list.


## License
This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.
If not, see http://www.gnu.org/licenses/.
