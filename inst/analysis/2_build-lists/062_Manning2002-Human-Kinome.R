#!/usr/bin/env Rscript

## ---- input 

root <- "~/projects/geneLists"
script <- "062_Manning2002-Human-Kinome" # no file-extension
short_name <- "hsKinome"
data_source <- "https://pubmed.ncbi.nlm.nih.gov/12471243/"
data_url <- "http://kinase.com/static/colt/data/human/kinome/tables/Kincat_Hsap.08.02.xls"

## The best article identifier is its DOI!
# * DOI: 10.1126/science.1075762

## ---- documentation

## FIXME: clean-up and shorten docs derived directly from Manning2002

# data are from 'Table S1', Manning2002 (Science):

###############################################################################
## Table S1: Sequences and Classification of Human Kinases
#   This table lists the classification and sequences for all human protein
# kinase genes and pseudogenes.  
#
###############################################################################
## Atypical Protein Kinases (aPKs)
###############################################################################
#   These proteins lack sequence similarity to the
# ePK domain HMM profile, but have been shown experimentally to have protein kinase
# activity, or are clear homologs of aPKs with demonstrated protein kinase
# activity. All aPK families are small, several having just one member
# invertebrates, and none in invertebrates. Other aPK families may remain to be
# discovered bybiochemical methods, but since most atypical families are small,
# and most biochemically-discovered kinases are ePKs, it is unlikely that many new
# atypical kinases will be discovered. Conversely, some of these aPKs may be false
# positives, if the reports of kinase activity are not correct. 
###############################################################################
## Alpha Kinases
###############################################################################
#   The progenitors of this family are the myosin heavy chain kinases (MHCKs) of
# Dictyosteliumdiscoideum. While these are evolutionarily restricted, they are
# similar to the eukaryoticelongation factor 2 kinase (eEF2K) found in most
# eukaryotes. Several other mammalian geneshave been found to be homologous to
# these, including the channel kinases Chak1 and Chak2,which are multi-pass
# transmembrane proteins which act as kinases and as ion channels.
###############################################################################
# Multiple
# sequence alignment shows that thePI34K domains of PIKKs form a distinct domain
# subfamily. The PI34K domain is structurallyrelated to the ePK domain; a
# structural alignment of ePK and PI34K domains shows similarstructure and
# conservation of most of the catalytically conserved residues, including
# the'catalytic' K, the DFG motif and the HRD motif (modified to DRH) 
# (see for instance Walker etal, 1999). A6This family consists of human A6 and A6r
# genes, along with homologs in fly, worm and
# yeast.
###############################################################################
# A6 was first cloned as a phospho-protein by Beeler et al (1994), who
# demonstrated tyrosinekinase activity by bacterial expression A6 fusion protein
# in an in vitro kinase activity. A6r wasdiscovered by Rohwer et al (1999) by
# two-hybrid interaction with PKC zeta; they show that bothA6 and A6r bind ATP,
# however, they failed to see kinase activity by either protein. To the bestof our
# knowledge, only one report supports the A6 family being a family of protein
# kinases.ABC1/ADCK (ABC1 domain containing kinase)This conserved family was
# identified as putative kinases by sequence alignment methods 
# (Psi-Blast and HMMs) which show a
# domain that is weakly similar to the ePK domain, withparticular conservation of
# the most conserved catalytic motifs. Their kinase similarity was firstpublished
# by Leonard et al in 1998. Despite the lack of overall sequence conservation with
# theePK domain, these kinases contain candidates for the most conserved kinase
# motifs, including the VAIK catalytic motif (VAVK, VAMK), the DFG motif, and a
# QTD motif that may take theplace of the HRD motif.PDK - Pyruvate Dehydrogenase
# KinasesThis family of mitochondrial kinases contains a domain which is similar
# to prokaryotic histidinekinases, but has been biochemically to phosphorylate
# serine rather than histidine. Crystalstructures (Machius et al, 2001, Steussy et
# 						 al, 2001) confirm that the PDK
# domain fold is similarto that of histidine kinases and Hsp90, and is distinct
# from the ePK domain.RIOThis family has 3 clear subfamilies, with one member of
# each in fly, worm and human. Yeast hastwo members (in the RIO1 and RIO2
# 						   subfamilies) and the fungus
# Aspergillus nidulans has amember of the third subfamily, RIO3. Homologs are also
# present in several archeal genomes.Yeast RIO1 was recently published to have
# serine kinase activity by Angermayr et al (2002). Thesequences do not align with
# the eukaryotic protein kinase domain, but many of the catalyticresidues are
# strongly conserved in the RIO family, and overall structural similarity to ePKs
# hasbeen predicted by Leonard et al (1998)BRD - Bromodomain KinasesThis family
# consists of the BRD2 (RING3) kinase and homologs in human and modelorganisms.
# Dennis and Green (1996) first identified BRD2 as an autophosphorylating
# nuclear-specific protein in Hela extracts. Recombinant BRD2 expressed in E. coli
# showed kinase activityin an in vitro assay, which was abolished by mutation of
# the putative catalytic lysine (K578A).The kinase activity was only seen when the
# purified recombinant protein was first incubated withHela cell extract and
# repurified, possibly due to the need for activation of BRD2 byphosphorylation by
# another kinase. Recombinant BRD2 purified from COS cells also showed invitro
# kinase activity, without the need for an activation step. Alignment of BRD2 with
# its humanhomologs (BRDT, BRD3, BRD4), Drosophila fsh and an anonymous worm
# homolog shows thepresence of two bromodomains, and a third conserved region,
# which contains some similarity tothe ePK domain by secondary structure
# prediction (I. Grigoriev, unpublished).TAF ñ TATA binding factor associated
# factorsTAF1 (TAF II-250) is a component of the basal transcriptional machinery,
# and exists as a singlecopy gene in all fully-sequenced eukaryotes. It has no
# close homologs. Dikstein et al (1996)report that TAF1 is a protein kinase and
# contains two regions which can independentlyphosphorylate the basal
# transcription factor RAP74. In vitro kinase assays were carried out
# withimmunopurified TFIID from Hela cells or with cloned TAF1 transfected into
# insect Sf9 cells orE. coli. Deletion mapping showed that two independent
# regions, each less than 470 AA long hadkinase activity, though neither had
# significant sequence similarity to each other, to proteinkinases, or to any
# other proteins. Later studies (Solow et al, 2001, OíBrien and Tijan,
# 			       1998)confirm the result and perform a finer
# mapping of the N-terminal kinase region. In 2002, Wangand Page revelated the
# presence of TAF1L, a retrotransposed copy of TAF1 present in humanand old-world
# primates, which is expressed during spermatogenesis and substitutes for TAF1 ina
# cellular assay.  BCRBest known as the fusion partner of the Abl kinase in
# chronic myologenous leukemia, the BCRgene itself also has protein kinase
# activity. Maru and Witte (1991) showed that highly-purifiedBCR has auto- and
# trans-phosphorylation activities, and later mapped cysteines and tyrosinesthat
# were critical for kinase activity. The kinase domain appears to be a recent
# addition to theprotein: the human ABR gene is about 70% identical in AA
# sequence, but lacks the N-terminalputative kinase domain. The Drosophila gene
# EG:23E12.2 (gi|7289304) also lacks the N-terminal putative kinase domain, but is
# 36% identical over 750 AA in the C-terminus. The regionconserved between these
# three proteins includes a GTPase activator domain.FASTKThe human Fas-activated
# s/t kinase (FASTK) was characterised by Tian et al (1995) as a kinasewhich was
# dephosphorylated and activated by Fas-mediated apoptosis. The nuclear TIA-1
# RNA-binding protein, a putative apoptosis effector, was identified as a
# substrate. Differentialexpression of FASTK in apoptotic cells has also been
# reported by Brutsche et al (2001). A closemouse ortholog is known (NP_148936.1),
# with 89% AA identity. Two other very distant putativemammalian homologs have
# been sequenced, but with ~27% identity between these genes, arenot close enough
# to confidently assign a kinase function to them.G11This family consists of a
# single gene called G11 or STK19, which was shown by Gomez-Escobaret al (1998).
# to have serine/threonine kinase activity against alpha casein and histone,
# showedATP-binding function and identified a putative catalytic lysine required
# for function (unlikeePKs, this lysine is near the C-terminal end of the
# 	      protein). The G11 protein was made intransfected insect or COS-7
# cells and immunopreciptated, so it remains possible that the kinaseactivity was
# due to a tightly bound protein.Clear homologs are found in rat and mouse 
# (~85% AA identity throughout),
# and a divergentputative homolog fragment has been sequenced in zebrafish (45%
# identity over 57 AA), but noobvious homolog has been seen in
# any other organism.
 
###############################################################################
## TIF1 Transcriptional Intermediary Factor 1 family
###############################################################################
#   A family of three Transcriptional Intermediary Factor 1 genes (TIF1a,b,g),
#   of which TIF1a has been shown to be a protein kinase (Fraser at al, 1998);
#   the two other genes are similar across their full length and so also likely
#   to be kinases.  Life the TAFs, TIFs are involved in thetranscriptional
#   machinery, and are though to phorphorylate several other
#   TATA-associatedfactors. Also similar to the TAFs and to BRD kinases, the
#   TIFs contain bromodomains,suggesting that they may be in some way inovlved
#   in the kinase functions of these proteins. 
###############################################################################
## H11
###############################################################################
# The H11 gene is a homolog of the ICP10 gene of Herpes simpex virus, both of
# which have been indicated to have kinase activity (Nelson et al, 1996, Smith et
# al, 2000).  H11 kinase activity wasdemonstrated in protein immnopurified from
# bacterial and eukaryotic expression systems. Thekinase activity was lost when
# a putative catlytic lysine was mutated. Smith et al noted weaksimilaritites
# between H11 and the FASTK atypical kinase, but it is not known if these
# aresignificant. H11 belongs to the HSP20/Alpha Crystallin family of proteins,
# and has a number ofclose and moderate homologs. Since the kinase region hasnít
# been mapped in this gene, it wasnot possible to confidently say if any of
# these homologs also have kinase activity, so only H11has been included in the
# kinome catalog.

## ---- prepare the renv

devtools::load_all(root, quiet=TRUE)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})


## ---- download the data

# downloaded on 01/20/2021
download.file(data_url,destfile=basename(data_url)) # 2.9 MB

# read with readxl
df <- readxl::read_xls(basename(data_url))

# rm file
invisible(file.remove(basename(data_url)))

# str(df), key columns:
# tibble [620 × 35] (S3: tbl_df/tbl/data.frame)
#   $ Name                                                   : chr [1:620] "AKT1" "AKT2" "AKT3" "CRIK" ...
#   $ SK                                                     : chr [1:620] "SK018" "SK019" "SK020" "SK695" ...
#   $ Group                                                  : chr [1:620] "AGC" "AGC" "AGC" "AGC" ...
#   $ Entrez_GeneID                                          : num [1:620] 207 208 10000 11113 1760 ...


## ---- load the data

entrez <- df$Entrez_GeneID
is_na <- is.na(entrez)
groups <- df$Group[!is_na]
gene_list <- split(entrez[!is_na],groups)

# drop na, get unique genes
gene_list <- lapply(gene_list, function(x) unique(x[!is.na(x)]))

# summarize
l <- sapply(gene_list,length) 
data.table(group=names(l),k=l) %>% knitr::kable()


## ---- save gene_list

gmtdir <- file.path(root,"gmt")
datadir <- file.path(root,"data")
funcdir <- file.path(root,"R")
myfile <- file.path(gmtdir, script)

write_gmt(gene_list,data_source, myfile)

documentDataset(myfile, short_name, funcdir,datadir)
