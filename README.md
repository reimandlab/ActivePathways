# ActivePathways

**July 10th 2020: ActivePathways is now available in CRAN: https://cran.r-project.org/web/packages/ActivePathways/index.html.**

ActivePathways is a tool for multivariate pathway enrichment analysis. Pathway enrichment analysis identifies gene sets, such as pathways or Gene Ontology terms, that are over-represented in a list of genes of interest. ActivePathways uses a data fusion method to combine multiple omics datasets, prioritizes genes based on the significance of signals from the omics datasets, and performs pathway enrichment analysis of these prioritized genes. Using this strategy, we can find pathways and genes supported by single or multiple omics datasets, as well as additional genes and pathways that are only apparent through data integration and remain undetected in any single dataset alone. 

ActivePathways is published in Nature Communications with the PCAWG Pan-Cancer project. 

Marta Paczkowska^, Jonathan Barenboim^, Nardnisa Sintupisut, Natalie S. Fox, Helen Zhu, Diala Abd-Rabbo, Miles W. Mee, Paul C. Boutros, PCAWG Drivers and Functional Interpretation Working Group, Jüri Reimand & PCAWG Consortium. Integrative pathway enrichment analysis of multivariate omics data. *Nature Communications* 11 735 (2020) (^ - co-first authors)
https://www.nature.com/articles/s41467-019-13983-9
https://www.ncbi.nlm.nih.gov/pubmed/32024846 

## Installation

#### From CRAN
Open R and run `install.packages('ActivePathways')`

#### Using devtools on our GitHub repository
Using the R package `devtools`, run
`devtools::install_github('https://github.com/reimandlab/ActivePathways', build_vignettes = TRUE)`

#### From source on our GitHub repository
Clone the repository, for example using `git clone https://github.com/reimandlab/ActivePathways.git`. 

Open R in the directory where you cloned the package and run `install.packages("ActivePathways", repos = NULL, type "source")`



## Using ActivePathways

* See the vignette for more details. Run `browseVignettes(package='ActivePathways')` in R.


### Examples
The simplest use of ActivePathways requires only a data table (matrix of p-values) and a list of gene sets in the form of a GMT [(Gene Matrix Transposed)](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29) file. The data table must be in the form of numerical matrix and cannot contain any missing values.
```R

library(ActivePathways)

##
# Run an example using the data files included in the ActivePathways package. 
##

fname_scores <- system.file("extdata", "Adenocarcinoma_scores_subset.tsv", package = "ActivePathways")
fname_GMT <- system.file("extdata", "hsapiens_REAC_subset.gmt", package = "ActivePathways")

##
# Numeric matrix of p-values is required as input. 
# NA values are converted to P = 1.
##

scores <- read.table(fname_scores, header = TRUE, row.names = 'Gene')
scores <- as.matrix(scores)
scores[is.na(scores)] <- 1


##
# Main call of ActivePathways function:
##

enriched_pathways <- ActivePathways(scores, fname_GMT) 

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
# Show enriched genes of the first pathway 'DAP12 signalling' 
# the column `overlap` displays genes of the integrated dataset (from 
# data fusion, i.e., p-value merging) that occur in the given pathway.
# Genes are ranked by joint significance across input omics datasets.
##

enriched_pathways[["overlap"]][[1]]
# [1] "TP53"   "PIK3CA" "KRAS"   "PTEN"   "BRAF"   "NRAS"   "B2M"    "CALM2"
# [9] "CDKN1A" "CDKN1B"

##
# Save the resulting pathways as a Comma-Separated Values (CSV) file for spreadsheets 
#  and computational pipelines.
# the data.table object cannot be saved directly as text.
##

export_as_CSV(enriched_pathways, "enriched_pathways.csv")


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


# Visualising pathway enrichment results using enrichment maps in Cytoscape

The Cytoscape software and the EnrichmentMap app provide powerful tools to visualise the enriched pathways from `ActivePathways` as a network (i.e., an Enrichment Map). To facilitate this visualisation step, `ActivePathways` provides the files needed for building enrichment maps. To create these files, a file prefix must be supplied to `ActivePathways` using the argument `cytoscape.file.tag`. The prefix can be a path to an existing writable directory.
 
```{r}
res <- ActivePathways(scores, gmt.file, cytoscape.file.tag = "enrichmentMap__")
```
Four files are written using the prefix:

* `enrichmentMap__pathways.txt` contains the table of significant terms (i.e. molecular pathways, biological processes, other gene sets) and the associated adjusted P-values. Note that only terms with `adjusted.p.val <= significant` are written.

* `enrichmentMap__subgroups.txt` contains a matrix indicating the columns of the input matrix of P-values that contributed to the discovery of the corresponding pathways. These values correspond to the `evidence` evaluation of input omics datasets discussed above, where a value of one indicates that the pathway was also detectable using a specific input omics dataset. A value of zero indicates otherwise. This file will be not generated if a single-column matrix of scores corresponding to just one omics dataset is provided to `ActivePathways`.

* `enrichmentMap__pathways.gmt` contains a shortened version of the supplied GMT file which consists of only the significant pathways detected by `ActivePathways`. 

* `enrichmentMap__legend.pdf` is a pdf file that displays a color legend of different omics datasets visualised in the enrichment map that can be used as a reference to the generated enrichment map.

## Creating enrichment maps using results of ActivePathways 

Pathway enrichment analysis often leads to complex and redundant results. Enrichment maps are network-based visualisations of pathway enrichment analyses. Enrichment maps can be generated in the Cytoscape software using the EnrichmentMap app. The enhancedGraphics app is also required. See the vignette for details: `browseVignettes(package='ActivePathways')`.


## Required software

1.	Cytoscape, see <https://cytoscape.org/download.html>
2.	EnrichmentMap app of Cytoscape, see menu Apps>App manager or <http://apps.cytoscape.org/apps/enrichmentmap> 
3.	EhancedGraphics app of Cytoscape, see menu Apps>App manager or <http://apps.cytoscape.org/apps/enhancedGraphics> 

## Creating the enrichment map

* Open the Cytoscape software. 
* Select *Apps -> EnrichmentMap*. 
* In the following dialogue, click the button `+` *Add Data Set from Files* in the top left corner of the dialogue.
* Change the Analysis Type to Generic/gProfiler/Enrichr.
* Upload the files `enrichmentMap__pathways.txt` and `enrichmentMap__pathways.gmt` in the *Enrichments* and *GMT* fields, respectively. 
* Click the checkbox *Show Advanced Options* and set *Cutoff* to 0.6.
* Then click *Build* in the bottom-right corner to create the enrichment map. 

![](https://github.com/reimandlab/ActivePathways/blob/master/vignettes/CreateEnrichmentMapDialogue_V2.png)

![](https://github.com/reimandlab/ActivePathways/blob/master/vignettes/NetworkStep1_V2.png)


## Colour the nodes of the network to visualise supporting omics datasets

To color nodes in the network (i.e., molecular pathways, biological processes) according to the omics datasets supporting the enrichments, the third file `enrichmentMap__subgroups.txt` needs to be imported to Cytoscape directly. To import the file, activate the menu option *File -> Import -> Table from File* and select the file `enrichmentMap__subgroups.txt`. In the following dialogue, select *To a Network Collection* in the dropdown menu *Where to Import Table data*. Click OK to proceed. 

![](https://github.com/reimandlab/ActivePathways/blob/master/vignettes/ImportStep_V2.png)

Next, Cytoscape needs to use the imported information to color nodes using a pie chart visualisation. To enable this click the Style tab in the left control panel and select the Image/Chart1 Property in a series of dropdown menus (*Properties -> Paint -> Custom Paint 1 -> Image/Chart 1*). 

![](https://github.com/reimandlab/ActivePathways/blob/master/vignettes/PropertiesDropDown2_V2.png)

The *image/Chart 1* property now appears in the Style control panel. Click the triangle on the right, then set the *Column* to *instruct* and the *Mapping type* to *Passthrough*. 

![](https://github.com/reimandlab/ActivePathways/blob/master/vignettes/StylePanel_V2.png)

This step colours the nodes corresponding to the enriched pathways according to the supporting omics datasets, based on in the scores matrix initially analysed in `ActivePathways`. 

![](https://github.com/reimandlab/ActivePathways/blob/master/vignettes/NetworkStep2_V2.png)

To allow better interpretation of the enrichment map, `ActivePathways` generates a color legend in the file `enrichmentMap__legend.pdf` that shows which colors correspond to which omics datasets. 

![](https://github.com/reimandlab/ActivePathways/blob/master/vignettes/LegendView.png)

Note that one of the colors corresponds to a subset of enriched pathways with *combined* evidence that were only detected through data fusion and P-value merging and not when any of the input datasets were detected separately. This exemplifies the added value of integrative multi-omics pathway enrichment analysis. 

## References

* See the vignette for more details: `browseVignettes(package='ActivePathways')`.

* Integrative Pathway Enrichment Analysis of Multivariate Omics Data. Paczkowska M, Barenboim J, Sintupisut N, Fox NS, Zhu H, Abd-Rabbo D, Mee MW, Boutros PC, PCAWG Drivers and Functional Interpretation Working Group; Reimand J, PCAWG Consortium. Nature Communications (2020) <https://pubmed.ncbi.nlm.nih.gov/32024846/> <https://doi.org/10.1038/s41467-019-13983-9>.

* Pathway Enrichment Analysis and Visualization of Omics Data Using g:Profiler, GSEA, Cytoscape and EnrichmentMap. Reimand J, Isserlin R, Voisin V, Kucera M, Tannus-Lopes C, Rostamianfar A, Wadi L, Meyer M, Wong J, Xu C, Merico D, Bader GD. Nature Protocols (2019) <https://pubmed.ncbi.nlm.nih.gov/30664679/> <https://doi.org/10.1038/s41596-018-0103-9>.



