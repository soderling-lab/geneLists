# Building Gene lists for enrichment analysis with the AnRichment Package.

## Datasets:
1. DBD - Developmental Brain Disorders (DBD) Genes Database
Description: A curated resource of genes associated with 6 DBDs (ID, autism, ADHD,
schizophrenia, bipolar disorder, and epilspy) compiled from the literature.
Dataset: "All"
Source: [url](http://dbd.geisingeradmi.org/#additional-information)
Mouse Genes: 443
Disorders: 6
Note:

2. DisGeneNET
Description: A database of gene-disease associations from public data sources
and the literature. Includes data from UniPrto, CGI, ClinGen, Genomics England
Panel App, PsyGeNET, Orphanet, HPO, and CTD (human). Genes compiled using the 
BeFree System are extracted from MEDLINE abstracts.
Dataset: "All Disease Genes : Mental or Behavioral Dysfunction/Mental Proces"
Source: [url](https://www.disgenet.org/static/disgenet_ap1/files/downloads)
Human Genes: ???
Mouse Genes: 3,118
Disorders:  12,491 
Note: Removed gene groups smaller than 3 or greater than 500.

3. SFARI
Description: an online database of human genes associated with autism.
Dataset: "SFARI"
Source: [url](https://www.sfari.org/resource/sfari-gene/)
Human Genes: 909
Mouse Genes: 859
Disorders: Autism spectrum disorders.
Note:

4. SynGO
Description: 
Dataset: ""
Source: [url]()
Human Genes: 
Mouse Genes:
Disorders: 
Note:

5. URMC DBDB - University of Rochester Medical Center (URMC) Developmental Brain Disorders Database (DBDB)
Description: A databse of DBD-associated genes for clinician-scientists. 
Dataset: "Associations"
Source: [url]("https://www.dbdb.urmc.rochester.edu/associations/list")
Human Genes: 600
Mouse Genes: 568
Disorders: 49
Note: Removed 39 gene groups with less than 3 genes.

6. Gene2Phenotype DBD
Description: G2P is a publicly-accessible online system designed to facilitate the
development, validation, curation and distribution of large-scale,
evidence-based datasets for use in diagnostic variant filtering. Each G2P entry
associates an allelic requirement and a mutational consequence at a defined
locus with a disease entity. A confidence level and evidence link are assigned
to each entry. The Developmental Disease (DD) database was curated by Helen V Firth 
and David FitzPatrick
Dataset: "DD:Brain/Cognition"
Source: [url]("https://www.ebi.ac.uk/gene2phenotype/downloads")
Human Genes: 1992
Mouse Genes: 1161
Disorders: 15 
Note:
