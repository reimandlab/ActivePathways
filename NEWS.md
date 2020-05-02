# activePathways 1.0.2

## Major Changes

* Changed logic of creating files for Enrichment Map: the user can provide 
the parameter "cytoscape.file.tag" for creating the required files. If the 
parameter is NA (default), no files are created. No directories are created. 
* removing command line parameter "return.all" as it was redundant with the 
parameter "significance".
* removing parameter "reanalyze" to simplify the package and leave the structuring 
of results up to the user.
* removing dependency of metap R package. As a result, only Fisher's and Brown's p-value 
merging options are available.


# activePathways 1.0.1

## Major Changes

* Rename to activePathways

# activeDriverPW 1.0.0

## Major Changes
* Add a geneset.filter option to activeDriverPW to filter the GMT by geneset size
* Add vignettes
* Export merge_p_values method