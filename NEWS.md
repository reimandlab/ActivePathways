# ActivePathways 1.1

## Major Changes

### ActivePathways 1.1
* Updated the gmt filtering step in ActivePathways.R. If a custom gene background vector is provided, then gmt gene sets share zero overlap
with the background are removed. This step also reduces the genes in the retained gmt gene sets to those present in the background.

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
