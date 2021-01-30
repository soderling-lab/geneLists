#' uniprot_map
#'
#' gene mapping data from the Universal Protein Resource (UniProt)
#' 
#' uniprot_map - a data.frame object containing the following mappings: 
#' 1. UniProtKB-AC
#' 2. UniProtKB-ID
#' 3. GeneID (EntrezGene)
#' 4. RefSeq
#' 5. GI
#' 6. PDB
#' 7. GO
#' 8. UniRef100
#' 9. UniRef90
#' 10. UniRef50
#' 11. UniParc
#' 12. PIR
#' 13. NCBI-taxon
#' 14. MIM
#' 15. UniGene
#' 16. PubMed
#' 17. EMBL
#' 18. EMBL-CDS
#' 19. Ensembl
#' 20. Ensembl_TRS
#' 21. Ensembl_PRO
#' 22. Additional PubMed
#' NOTE: columns 'Accession' (UniProt Accession), 'ID' (Uniprot ID),
#' and 'Entrez' (GeneID), have been addded for ease-of-use and clarity.
#' 
#' @docType data
#'
#' @usage data()
#'
#' @format An object of class \code{"data.table"};
#'
#' @keywords datasets
#'
#' @references none
#'
#' @source \href{source}{ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/}

#' The Universal Protein Resource (UniProt), a collaboration between the European
#' Bioinformatics Institute (EBI), the SIB Swiss Institute of Bioinformatics, and
#' the Protein Information Resource (PIR), is comprised of three databases, each
#' optimized for different uses. The UniProt Knowledgebase (UniProtKB) is the
#' central access point for extensively curated protein information, including
#' function, classification and cross-references. The UniProt Reference Clusters
#' (UniRef) combine closely related sequences into a single record to speed up
#' sequence similarity searches. The UniProt Archive (UniParc) is a comprehensive
#' repository of all protein sequences, consisting only of unique identifiers and
#' sequences.
#'
#'   UniProt LICENSE
#' We have chosen to apply the Creative Commons Attribution (CC BY 4.0) License
#' (https://creativecommons.org/licenses/by/4.0/) to all copyrightable parts of
#' our databases.
#' 
#' (c) 2002-2020 UniProt Consortium
#' 
#'   UniProt DISCLAIMER
#' We make no warranties regarding the correctness of the data, and disclaim
#' liability for damages resulting from its use. We cannot provide unrestricted
#' permission regarding the use of the data, as some data may be covered by patents
#' or other rights.
#' 
#' Any medical or genetic information is provided for research, educational and
#' informational purposes only. It is not in any way intended to be used as a
#' substitute for professional medical advice, diagnosis, treatment or care.
