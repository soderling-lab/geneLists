https://www.syngoportal.org 
SynGO dataset version/release: 20180731

syngo_genes.xlsx contains all genes in the SynGO database.

syngo_annotations.xlsx lists all individual annotation statements.

syngo_ontologies.xlsx lists all ontology terms and their respective annotated genes.

In case you are matching/integrating this data with other GO resources, note that some SynGO ontology terms are not part of the GO database, their identifiers are easily recognised by the "SYNGO:" prefix as opposed to "GO:". Annotations against these terms is integrated into GO, but translated into an annotation+extensions as some of the SynGO terms do not comply with GO guidelines (eg; a term like "presynaptic processes" is too generic for GO as they'd end up with thousands of suchs terms).

The .json files are intended for bioinformatic application of the SynGO dataset. These contain the 'gene sets' for each SynGO ontology term. Multiple ID mapping types are available. Note that the 'direct' objects contain only annotations specifically against some term, while 'aggregate' contains direct annotations + child terms. For most bioinformatic applications, such as geneset enrichment analysis, you probably want to use the latter.
Do consider using the SynGO color-coding tool to visualize your own SynGO ontology term summary scores.