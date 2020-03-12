# activePathways

activePathways is a tool for multivariate pathway enrichment analysis. Pathway enrichment analysis identifies gene sets, such as pathways or Gene Ontology terms, that are over-represented in a list of genes of interest. activePathways implements an algorithm to perform pathway enrichment analysis on a data set that contains multiple variables for each gene. For example, the data may contain the confidence that the gene is a driver as reported by several tools, or it may contain differential expression values across different genetic regions such as coding regions and promoter regions.

ActivePathways is published in Nature Communications with the PCAWG Pan-Cancer project. 

Marta Paczkowska^, Jonathan Barenboim^, Nardnisa Sintupisut, Natalie S. Fox, Helen Zhu, Diala Abd-Rabbo, Miles W. Mee, Paul C. Boutros, PCAWG Drivers and Functional Interpretation Working Group, Jüri Reimand & PCAWG Consortium. Integrative pathway enrichment analysis of multivariate omics data. *Nature Communications* 11 735 (2020) https://www.nature.com/articles/s41467-019-13983-9 https://www.ncbi.nlm.nih.gov/pubmed/32024846 (^ - co-first authors)

## Installation

#### devtools:
Using the R package `devtools`, run
`devtools::install_github('https://github.com/reimandlab/activePathways')`

#### From source
Clone the repository: `https://github.com/reimandlab/activePathways.git`
Open R in the directory you cloned the package in and run `install.packages('activePathways', repos=NULL)`

## Using activePathways

### Examples
The simplest use of activeDriver requires only a data table and a list of gene sets in the form of a GMT [(Gene Matrix Transposed)](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29) file. The data must be in the form of numerical matrix and cannot contain any missing values
```
scores <- read.table('example_data.txt', header=TRUE, row.names='Gene')
scores <- as.matrix(scores)
scores[is.na(scores)] <- 1
scores
#                      X3UTR        X5UTR          CDS     promCore
## A2M           1.0000000000 6.679353e-01 9.051708e-01 4.499201e-01
## AAAS          1.0000000000 8.501202e-01 7.047723e-01 7.257641e-01
## ABAT          0.9664125635 8.405470e-02 7.600985e-01 1.903789e-01
## ABCC1         0.9383431571 9.198887e-01 2.599319e-01 2.980455e-01
##  [ reached getOption("max.print") -- omitted 2410 rows ]

res <- activePathways(scores, 'example_genesets.gmt')
```

More thorough documentation of the activePathways function can be found in R with `?activePathways`, and complete tutorials can be found with `browseVignettes(package='activePathways')`
