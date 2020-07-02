# ActivePathways

**May 2nd 2020: ActivePathways is currently not available from CRAN due to broken dependencies. This GitHub repository provides quick fixes to a number of recently-emerged issues. An updated CRAN package (ActivePathways 1.0.2) will be available shortly. We apologize for the inconvenience and thank you for your patience while we fix the code and documentation.**

ActivePathways is a tool for multivariate pathway enrichment analysis. Pathway enrichment analysis identifies gene sets, such as pathways or Gene Ontology terms, that are over-represented in a list of genes of interest. ActivePathways uses a data fusion method to combine multiple omics datasets, prioritizes genes based on the significance of signals from the omics datasets, and performs pathway enrichment analysis of these prioritized genes. Using this strategy, we can find pathways and genes supported by single or multiple omics datasets, as well as additional genes and pathways that are only apparent through data integration and remain undetected in any single dataset alone. 

ActivePathways is published in Nature Communications with the PCAWG Pan-Cancer project. 

Marta Paczkowska^, Jonathan Barenboim^, Nardnisa Sintupisut, Natalie S. Fox, Helen Zhu, Diala Abd-Rabbo, Miles W. Mee, Paul C. Boutros, PCAWG Drivers and Functional Interpretation Working Group, Jüri Reimand & PCAWG Consortium. Integrative pathway enrichment analysis of multivariate omics data. *Nature Communications* 11 735 (2020) (^ - co-first authors)
https://www.nature.com/articles/s41467-019-13983-9
https://www.ncbi.nlm.nih.gov/pubmed/32024846 

## Installation

#### From CRAN (not available temporarily)
Open R and run `install.packages('ActivePathways')`

#### Using devtools on our GitHub repository
Using the R package `devtools`, run
`devtools::install_github('https://github.com/reimandlab/ActivePathways')`

#### From source on our GitHub repository
Clone the repository, for example using `git clone https://github.com/reimandlab/ActivePathways.git`. 

Open R in the directory where you cloned the package and run `install.packages("ActivePathways", repos = NULL, type "source")`



## Using ActivePathways

### Examples
The simplest use of ActivePathways requires only a data table (matrix of p-values) and a list of gene sets in the form of a GMT [(Gene Matrix Transposed)](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29) file. The data table must be in the form of numerical matrix and cannot contain any missing values.
```

library(ActivePathways)

##
# Run an example using the data files included in the ActivePathways package. 
##

fname_scores = system.file("extdata", "Adenocarcinoma_scores_subset.tsv", package = "ActivePathways")
fname_GMT = system.file("extdata", "hsapiens_REAC_subset.gmt", package = "ActivePathways")

##
# Numeric matrix of p-values is required as input. 
# NA values are converted to P = 1.
##

scores = read.table(fname_scores, header = TRUE, row.names = 'Gene')
scores <- as.matrix(scores)
scores[is.na(scores)] <- 1


##
# Main call of ActivePathways function:
##

enriched_pathways =  ActivePathways(scores, fname_GMT) 

#35 terms were removed from gmt because they did not make the geneset.filter
#91 rows were removed from scores because they are not found in the background


##
# list a few first results of the ActivePathways analysis
##

enriched_pathways[1:3,]

#        term.id         term.name adjusted.p.val term.size
#1: REAC:2424491   DAP12 signaling   4.491268e-05       358
#2:  REAC:422475     Axon guidance   2.028966e-02       555
#3:  REAC:177929 Signaling by EGFR   6.245734e-04       366
#                                   overlap       evidence
#1:     TP53,PIK3CA,KRAS,PTEN,BRAF,NRAS,...            CDS
#2: PIK3CA,KRAS,BRAF,NRAS,CALM2,RPS6KA3,... X3UTR,promCore
#3:     TP53,PIK3CA,KRAS,PTEN,BRAF,NRAS,...            CDS
#                            Genes_X3UTR Genes_X5UTR
#1:                                   NA          NA
#2: CALM2,ARPC2,RHOA,NUMB,CALM1,ACTB,...          NA
#3:                                   NA          NA
#                             Genes_CDS
#1: TP53,PTEN,KRAS,PIK3CA,BRAF,NRAS,...
#2:                                  NA
#3: TP53,PTEN,KRAS,PIK3CA,BRAF,NRAS,...
#                                Genes_promCore
#1:                                          NA
#2: EFNA1,IQGAP1,COL4A1,SCN2B,RPS6KA3,CALM2,...
#3:                                          NA


## 
# Examine a few lines of the two major types of input
##

##
# The scores matrix includes p-values for genes (rows) 
#   and evidence of different omics datasets (columns).
# This dataset includes predicted cancer driver mutations
#   in gene CDS, UTR and core promoter sequences
##

head(scores, n = 3)

#         X3UTR      X5UTR       CDS  promCore
#A2M  1.0000000 0.33396764 0.9051708 0.4499201
#AAAS 1.0000000 0.42506012 0.7047723 0.7257641
#ABAT 0.9664126 0.04202735 0.7600985 0.1903789

##
# GMT files include functional gene sets (pathways, processes).
# Each tab-separated line represents a gene set: 
#   gene set ID, description followed by gene symbols.
# Gene symbols in the scores table and the GMT file need to match. 
##

readLines(fname_GMT)[11:13]
#[1] "REAC:3656535\tTGFBR1 LBD Mutants in Cancer\tTGFB1\tFKBP1A\tTGFBR2\tTGFBR1\t"
#[2] "REAC:73927\tDepurination\tOGG1\tMPG\tMUTYH\t"
#[3] "REAC:5602410\tTLR3 deficiency - HSE\tTLR3\t" 



```

More thorough documentation of the ActivePathways function can be found in R with `?ActivePathways`, and complete tutorials can be found with `browseVignettes(package='ActivePathways')`.

### Creating enrichment maps in Cytoscape

In R, type in `browseVignettes(package='ActivePathways')` for a link to a tutorial with screenshots. 



