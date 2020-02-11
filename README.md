# gene-lists 
Gene lists for gene enrichment analysis.

## Developmental Brain Disorder (DBD)-Associated Genes
Human and animal genes associated with developmental brain disorders __(DBDs)__ or 
DBD phenotypes were compiled from several databases. 
Human genes were mapped to mouse homologs using the NCBI Homologene database.

With the exeption of SFARI gene groups, Disease/phenotype groups were 
filtered to remove groups with less than 3, or more than 500 genes.

In all, __1,306__ DBD-gene associations were compiled from __8__ databases
cooresponding to 984 unique DBD-associated mouse genes.

### 1. DBD Genes (Geisinger.org)
__Description:__ A curated resource of genes associated with 6 DBDs (ID, autism, ADHD,
schizophrenia, bipolar disorder, and epilspy) compiled from the literature.  
__Dataset:__ `'All'` -- includes LOF and missense variants.  
__Source:__ [url](http://dbd.geisingeradmi.org/#additional-information)  
__Mouse Genes:__ 443  
__Disorders/Phenotypes:__ 6  

### 2. DisGeneNET
__Description:__ A database of gene-disease associations from public data sources
and the literature. Includes data from UniProt, CGI, ClinGen, Genomics England
Panel App, PsyGeNET, Orphanet, HPO, and CTD. Genes compiled using the 
BeFree System are extracted from MEDLINE abstracts.  
__Dataset:__ `'All Disease Genes':'Mental or Behavioral Dysfunction'/'Mental Process'`  
__Source:__ [url](https://www.disgenet.org/static/disgenet_ap1/files/downloads)  
__Mouse Genes:__ 1718  
__Disorders/Phenotypes:__ 371  

### 3. SFARI Human Genes
__Description:__ an online database of human and animal genes associated with autism.  
__Dataset:__ `'SFARI'`  
__Source:__ [url](https://www.sfari.org/resource/sfari-gene/)  
__Mouse Genes:__ 859  
__Disorders/Phenotypes:__ 1  

### 4. SFARI Animal Genes
__Description:__ an online database of human and animal genes associated with autism.  
__Dataset:__ `'Animal'`  
__Source:__ [url](https://www.sfari.org/resource/sfari-gene/)  
__Mouse Genes:__ 189  
__Disorders/Phenotypes:__ 1  

### 5. URMC DBDB 
__Description:__ The University of Rochester Medical Center (URMC) 
Developmental Brain Disorders Database (DBDB) is a databse of 
DBD-associated genes intended for use by clinician-scientists.   
__Dataset:__ `'Associations'`  
__Source:__ [url]("https://www.dbdb.urmc.rochester.edu/associations/list")  
__Mouse Genes:__ 568  
__Disorders/Phenotypes:__ 49  

### 6. Gene2Phenotype DBD
__Description:__ G2P is an online database designed to facilitate the
development, validation, curation and distribution of large-scale,
evidence-based datasets for use in diagnostic variant filtering. Each G2P entry
associates an allelic requirement and a mutational consequence at a defined
locus with a disease entity. A confidence level and evidence link are assigned
to each entry.   
__Dataset:__ `'DD':'Brain/Cognition'`  
__Source:__ [url]("https://www.ebi.ac.uk/gene2phenotype/downloads")  
__Mouse Genes:__ 169  
__Disorders/Phenotypes:__ 15   
__Note:__ Genes associated with brain/congition phenotypes were kept.
The Developmental Disease (DD) database was curated by Helen 
V Firth and David FitzPatrick.   

### 7. Sanders et al., 2015.
__Description:__ 
__Dataset:__ `65genes_tadaFdrAscSscExomeSscAgpSmallDel +
59genes_tadaFdrAscSscExome`  
__Source:__ [url]("https://www.ncbi.nlm.nih.gov/pubmed/26402605")  
__Mouse Genes:__ 67  
__Disorders/Phenotypes:__ 1   
__Note:__ 

### 8. Satterstrom et al., 2015.
__Description:__   
__Dataset:__ `ASD/DDID`  
__Source:__ [url]("")  
__Mouse Genes:__ 96  
__Disorders/Phenotypes:__ 2   
__Note:__ 
