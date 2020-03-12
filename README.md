# ActivePathways

ActivePathways is a tool for multivariate pathway enrichment analysis. Pathway enrichment analysis identifies gene sets, such as pathways or Gene Ontology terms, that are over-represented in a list of genes of interest. ActivePathways uses a data fusion method to combine multiple omics datasets, prioritizes genes based on the significance of signals from the omics datasets, and performs pathway enrichment analysis of these prioritized genes. Using this strategy, we can find pathways and genes supported by single or multiple omics datasets, as well as additional genes and pathways that are only apparent through data integration and remain undetected in any single dataset alone. 

ActivePathways is published in Nature Communications with the PCAWG Pan-Cancer project. 

Marta Paczkowska^, Jonathan Barenboim^, Nardnisa Sintupisut, Natalie S. Fox, Helen Zhu, Diala Abd-Rabbo, Miles W. Mee, Paul C. Boutros, PCAWG Drivers and Functional Interpretation Working Group, Jüri Reimand & PCAWG Consortium. Integrative pathway enrichment analysis of multivariate omics data. *Nature Communications* 11 735 (2020) https://www.nature.com/articles/s41467-019-13983-9 https://www.ncbi.nlm.nih.gov/pubmed/32024846 (^ - co-first authors)

## Installation

#### From CRAN
Open R and run `install.packages('ActivePathways')`

#### Using devtools on our GitHub repository
Using the R package `devtools`, run
`devtools::install_github('https://github.com/reimandlab/activePathways')`

#### From source on our GitHub repository
Clone the repository: `https://github.com/reimandlab/activePathways.git`
Open R in the directory you cloned the package in and run `install.packages('activePathways', repos=NULL)`



## Using ActivePathways

### Examples
The simplest use of ActivePathways requires only a data table (matrix of p-values) and a list of gene sets in the form of a GMT [(Gene Matrix Transposed)](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29) file. The data table must be in the form of numerical matrix and cannot contain any missing values.
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

res <- ActivePathways(scores, 'example_genesets.gmt')
```

More thorough documentation of the ActivePathways function can be found in R with `?ActivePathways`, and complete tutorials can be found with `browseVignettes(package='ActivePathways')`.

### Creating enrichment maps in Cytoscape

In R, type in `browseVignettes(package='ActivePathways')` for a link to a tutorial with screenshots. 



