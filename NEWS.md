# ActivePathways 2.0.0

## Major Changes

### ActivePathways 2.0.1
* Changed how very small P-values are processed before P-value merging. P-values of '0' or anything less than '1e-300' are converted to '1e-300'.

### ActivePathways 2.0.0
* Incorporated scores_direction' and 'constraints_vector' parameters to ActivePathways() and merge_p_values() to account for the direction between datasets when performing p-value merging.
* Added the 'Stouffer' and 'Strube' p-value merging methods as alternatives to 'Fisher' and "Brown', respectively. 
* Changed the naming convention of parameters and objects, substituting a period '.' for an underscore '_'. 

### ActivePathways 1.1.1
* Added an option for the colours specified in the "custom_colors" parameter to be provided in any order, as long as the vector "names()" match the column names of the "scores" matrix. 
* Fixed an error where the "combined" contribution label was absent from the legend.pdf ActivePathways output file. 

### ActivePathways 1.1.0
* Updated the filtering procedure of gene sets in the GMT file when a custom gene background is provided. Given a background gene list, the GMT gene sets are first modified to only include the genes of the background list, and second, the gene sets are filtered by gene set size. Gene sets lacking any genes from the background list are removed. This update will result in a more lenient multiple testing correction in analyses with a custom background gene list.

### ActivePathways 1.0.4
* Added three new parameters to ActivePathways and prepareCystoscape functions. These include "color_palette", "custom_colors" and "color_integrated_only" to provide more options for node coloring in Cytoscape. 

### ActivePathways 1.0.3
* Removed dependency used for testing for CRAN compliance.

### ActivePathways 1.0.2
* Renamed package to ActivePathways from activePathways for consistency 
with function and publication
* Added new function export_as_CSV(res, file_name) to save data in 
spreadsheet-friendly formats
* Updated README-file with an actionable step-by-step tutorial
* Changed logic of creating files for Enrichment Map: the user can provide 
the parameter "cytoscape.file.tag" for creating the required files. If the 
parameter is NA (default), no files are created. No directories are created. 
* Removed the parameter "return.all" as it was redundant with the 
parameter "significance".
* Removed the parameter "reanalyze" to simplify the package and leave the structuring 
of results up to the user.
* Removed the dependency on the R package metap. As a result, only Fisher's and Brown's p-value 
merging options are available.
* Updated the vignette that now describes the ActivePathways package as well as the following steps of visualising results as enrichment maps in Cytoscape. 
