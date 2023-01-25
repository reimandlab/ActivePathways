#' ActivePathways
#'
#' @param scores A numerical matrix of p-values where each row is a gene and
#'   each column represents an omics dataset (evidence). Rownames correspond to the genes 
#'   and colnames to the datasets. All values must be 0<=p<=1. We recommend converting 
#'   missing values to ones. 
#' @param gmt A GMT object to be used for enrichment analysis. If a filename, a
#'   GMT object will be read from the file.
#' @param background A character vector of gene names to be used as a
#'   statistical background. By default, the background is all genes that appear
#'   in \code{gmt}.
#' @param geneset_filter A numeric vector of length two giving the lower and 
#'   upper limits for the size of the annotated geneset to pathways in gmt.
#'   Pathways with a geneset shorter than \code{geneset_filter[1]} or longer
#'   than \code{geneset_filter[2]} will be removed. Set either value to NA
#'   to not enforce a minimum or maximum value, or set \code{geneset_filter} to 
#'   \code{NULL} to skip filtering.
#' @param cutoff A maximum merged p-value for a gene to be used for analysis.
#'   Any genes with merged, unadjusted \code{p > significant} will be discarded 
#'   before testing.
#' @param significant Significance cutoff for selecting enriched pathways. Pathways with
#'   \code{adjusted_p_val <= significant} will be selected as results.
#' @param merge_method Statistical method to merge p-values. See section on Merging P-Values
#' @param correction_method Statistical method to correct p-values. See
#'   \code{\link[stats]{p.adjust}} for details.
#' @param cytoscape_file_tag The directory and/or file prefix to which the output files
#'   for generating enrichment maps should be written. If NA, files will not be written. 
#' @param color_palette Color palette from RColorBrewer::brewer.pal to color each
#'   column in the scores matrix. If NULL grDevices::rainbow is used by default.
#' @param custom_colors A character vector of custom colors for each column in the scores matrix.
#' @param color_integrated_only A character vector of length 1 specifying the color of the 
#'   "combined" pathway contribution.
#' @param scores_direction A numerical matrix of log2 transformed fold-change values where each row is a
#'   gene and each column represents a dataset (evidence). Rownames correspond to the genes
#'   and colnames to the datasets. We recommend converting missing values to ones. 
#'   Must contain the same dimensions as the scores parameter. Datasets without directional information should be set to 0.
#' @param expected_direction A numerical vector of +1 or -1 values corresponding to the expected
#'   directional relationship between columns in scores_direction. Datasets without directional information should
#'   be set to 0.
#'
#' @return A data.table of terms (enriched pathways) containing the following columns:
#'   \describe{
#'     \item{term_id}{The database ID of the term}
#'     \item{term_name}{The full name of the term}
#'     \item{adjusted_p_val}{The associated p-value, adjusted for multiple testing}
#'     \item{term_size}{The number of genes annotated to the term}
#'     \item{overlap}{A character vector of the genes enriched in the term}
#'     \item{evidence}{Columns of \code{scores} (i.e., omics datasets) that contributed 
#'          individually to the enrichment of the term. Each input column is evaluated 
#'          separately for enrichments and added to the evidence if the term is found.}
#'   }
#'
#' @section Merging P-values:
#' To obtain a single p-value for each gene across the multiple omics datasets considered, 
#' the p-values in \code{scores} #' are merged row-wise using a data fusion approach of p-value merging. 
#' The four available methods are:
#' \describe{
#'  \item{Fisher}{Fisher's method assumes p-values are uniformly
#'  distributed and performs a chi-squared test on the statistic sum(-2 log(p)).
#'  This method is most appropriate when the columns in \code{scores} are
#'  independent.}
#'  \item{Brown}{Brown's method extends Fisher's method by accounting for the
#'  covariance in the columns of \code{scores}. It is more appropriate when the
#'  tests of significance used to create the columns in \code{scores} are not
#'  necessarily independent. The Brown's method is therefore recommended for 
#'  many omics integration approaches.}
#'  \item{Stouffer}{Stouffer's method assumes p-values are uniformly distributed
#'  and transforms p-values into a Z-score using the cumulative distribution function of a
#'  standard normal distribution. This method is appropriate when the columns in \code{scores}
#'   are independent.}
#'  \item{Strube}{Strube's method extends Stouffer's method by accounting for the 
#'  covariance in the columns of \code{scores}.}
#' }
#'
#' @section Cytoscape:
#'   To visualize and interpret enriched pathways, ActivePathways provides an option
#'   to further analyse results as enrichment maps in the Cytoscape software. 
#'   If \code{!is.na(cytoscape_file_tag)}, four files will be written that can be used 
#'   to build enrichment maps. This requires the EnrichmentMap and enhancedGraphics apps.
#'
#' The four files written are:
#'   \describe{
#'     \item{pathways.txt}{A list of significant terms and the
#'     associated p-value. Only terms with \code{adjusted_p_val <= significant} are
#'     written to this file.}
#'     \item{subgroups.txt}{A matrix indicating whether the significant terms (pathways)
#'     were also found to be significant when considering only one column from
#'     \code{scores}. A one indicates that term was found to be significant 
#' 			when only p-values in that column were used to select genes.}
#'     \item{pathways.gmt}{A Shortened version of the supplied GMT
#'     file, containing only the significantly enriched terms in pathways.txt. }
#'     \item{legend.pdf}{A legend with colours matching contributions
#'     from columns in \code{scores}.}
#'   }
#'
#'   How to use: Create an enrichment map in Cytoscape with the file of terms
#'   (pathways.txt) and the shortened gmt file
#'   (pathways.gmt). Upload the subgroups file (subgroups.txt) as a table
#'   using the menu File > Import > Table from File. To paint nodes according 
#'   to the type of supporting evidence, use the 'style'
#'   panel, set image/Chart1 to use the column `instruct` and the passthrough
#'   mapping type. Make sure the app enhancedGraphics is installed. 
#'   Lastly, use the file legend.pdf as a reference for colors in the enrichment map.
#'
#' @examples
#'     fname_scores <- system.file("extdata", "Adenocarcinoma_scores_subset.tsv", 
#'          package = "ActivePathways")
#'     fname_GMT = system.file("extdata", "hsapiens_REAC_subset.gmt",
#'          package = "ActivePathways")
#'
#'     dat <- as.matrix(read.table(fname_scores, header = TRUE, row.names = 'Gene'))
#'     dat[is.na(dat)] <- 1
#'
#'     ActivePathways(dat, fname_GMT)
#'
#' @import data.table
#'
#' @export

ActivePathways <-  function(scores, gmt, background = makeBackground(gmt),
                            geneset_filter = c(5, 1000), cutoff = 0.1, significant = 0.05,
                            merge_method = c("Brown", "Fisher", "Stouffer","Strube"),
                            correction_method = c("holm", "fdr", "hochberg", "hommel",
                                                  "bonferroni", "BH", "BY", "none"),
                            cytoscape_file_tag = NA, color_palette = NULL, custom_colors = NULL, 
                            color_integrated_only = "#FFFFF0", scores_direction = NULL, 
                            expected_direction = NULL) {
  
  merge_method <- match.arg(merge_method)
  correction_method <- match.arg(correction_method)
  
  ##### Validation #####
  # scores
  if (!(is.matrix(scores) && is.numeric(scores))) stop("scores must be a numeric matrix")
  if (any(is.na(scores))) stop("scores may not contain missing values")
  if (any(scores < 0) || any(scores > 1)) stop("All values in scores must be in [0,1]")
  if (any(duplicated(rownames(scores)))) stop("scores matrix contains duplicated genes - rownames must be unique")
  
  # scores_direction and expected_direction
  if (xor(!is.null(scores_direction),!is.null(expected_direction))) stop("Both scores_direction and expected_direction must be provided")
  if (!is.null(scores_direction) && !is.null(expected_direction)){
        if (!(is.numeric(expected_direction) && is.vector(expected_direction))) stop("expected_direction must be a numeric vector")
        if (any(!expected_direction %in% c(1,-1,0))) stop("expected_direction must contain the values: 1, -1 or 0")
        if (!(is.matrix(scores_direction) && is.numeric(scores_direction))) stop("scores_direction must be a numeric matrix")
        if (any(is.na(scores_direction))) stop("scores_direction may not contain missing values")
        if (any(!rownames(scores_direction) %in% rownames(scores))) stop ("scores_direction gene names must match scores genes")
        if (is.null(colnames(scores)) || is.null(colnames(scores_direction))) stop("column names must be provided to scores and scores_direction")
        if (any(!colnames(scores_direction) %in% colnames(scores))) stop("scores_direction column names must match scores column names")
        if (length(expected_direction) != length(colnames(scores_direction))) stop("expected_direction should have the same number of entries as columns in scores_direction")
        if (any(expected_direction %in% 0) &&  !all(scores_direction[,expected_direction %in% 0] == 0)) 
              stop("scores_direction entries must be set to 0's for columns that do not contain directional information")
        if (!is.null(names(expected_direction))){
              if (!all.equal(names(expected_direction), colnames(scores_direction), colnames(scores)) == TRUE){
                    stop("the expected_direction entries should match the order of scores and scores_direction columns")
              }}}
        
  # cutoff and significant
  stopifnot(length(cutoff) == 1)
  stopifnot(is.numeric(cutoff))
  if (cutoff < 0 || cutoff > 1) stop("cutoff must be a value in [0,1]")
  stopifnot(length(significant) == 1)
  stopifnot(is.numeric(significant))
  if (significant < 0 || significant > 1) stop("significant must be a value in [0,1]")
  
  # gmt
  if (!is.GMT(gmt)) gmt <- read.GMT(gmt)
  if (length(gmt) == 0) stop("No pathways in gmt made the geneset_filter")
  if (!(is.character(background) && is.vector(background))) {
    stop("background must be a character vector")
  } 
  
  # geneset_filter
  if (!is.null(geneset_filter)) {
    if (!(is.numeric(geneset_filter) && is.vector(geneset_filter))) {
      stop("geneset_filter must be a numeric vector")
    }
    if (length(geneset_filter) != 2) stop("geneset_filter must be length 2")
    if (!is.numeric(geneset_filter)) stop("geneset_filter must be numeric")
    if (any(geneset_filter < 0, na.rm=TRUE)) stop("geneset_filter limits must be positive")
  }
  
  # custom_colors
  if (!is.null(custom_colors)){
    if(!(is.character(custom_colors) && is.vector(custom_colors))){
      stop("colors must be provided as a character vector")   
    } 
    if(length(colnames(scores)) != length(custom_colors)) stop("incorrect number of colors is provided")
  }
  if (!is.null(custom_colors) & !is.null(color_palette)){
    stop("Both custom_colors and color_palette are provided. Specify only one of these parameters for node coloring.")
  }
  
  if (!is.null(names(custom_colors))){
        if (!all(names(custom_colors) %in% colnames(scores))){
              stop("names() of the custom colors vector should match the scores column names")
        }
  }
  
  # color_palette
  if (!is.null(color_palette)){
    if (!(color_palette %in% rownames(RColorBrewer::brewer.pal.info))) stop("palette must be from the RColorBrewer package")
  }
  
  # color_integrated_only
  if(!(is.character(color_integrated_only) && is.vector(color_integrated_only))){
      stop("color must be provided as a character vector")   
  } 
  if(1 != length(color_integrated_only)) stop("only a single color must be specified")
	
  # contribution
  contribution <- TRUE
  if (ncol(scores) == 1) {
    contribution <- FALSE
    message("scores matrix contains only one column. Column contributions will not be calculated")
  }
  
  ##### filtering and sorting ####
    
  # Remove any genes not found in the background
  orig_length <- nrow(scores)
  scores <- scores[rownames(scores) %in% background, , drop=FALSE]
  if(!is.null(scores_direction)){
        scores_direction <- scores_direction[rownames(scores_direction) %in% background, , drop=FALSE]
  }
  if (nrow(scores) == 0) {
    stop("scores does not contain any genes in the background")
  }
  if (nrow(scores) < orig_length) {
    message(paste(orig_length - nrow(scores), "rows were removed from scores",
                  "because they are not found in the background"))
  }
  
	
  # Filter the GMT
  if (!all(background %in% unique(unlist(sapply(gmt, "[", c(3)))))){
        background_genes <- lapply(sapply(gmt, "[", c(3)), intersect, background)
        background_genes <- background_genes[lapply(background_genes,length) > 0]
        gmt <- gmt[names(sapply(gmt,"[",c(3))) %in% names(background_genes)]
        for (i in 1:length(gmt)) {
          gmt[[i]]$genes <- background_genes[[i]]
        }
  }
	
  if(!is.null(geneset_filter)) {
    orig_length <- length(gmt)
    if (!is.na(geneset_filter[1])) {
      gmt <- Filter(function(x) length(x$genes) >= geneset_filter[1], gmt)
    }
    if (!is.na(geneset_filter[2])) {
      gmt <- Filter(function(x) length(x$genes) <= geneset_filter[2], gmt)
    }
    if (length(gmt) == 0) stop("No pathways in gmt made the geneset_filter")
    if (length(gmt) < orig_length) {
      message(paste(orig_length - length(gmt), "terms were removed from gmt", 
                    "because they did not make the geneset_filter"))
    }
  }
  
  # merge p-values to get a single score for each gene and remove any genes
  # that don't make the cutoff
  merged_scores <- merge_p_values(scores, merge_method,scores_direction,expected_direction)
  merged_scores <- merged_scores[merged_scores <= cutoff]
  
  if (length(merged_scores) == 0) stop("No genes made the cutoff")
  
  # Sort genes by p-value
  ordered_scores <- names(merged_scores)[order(merged_scores)]
  
  ##### enrichmentAnalysis and column contribution #####
  
  res <- enrichmentAnalysis(ordered_scores, gmt, background)
  adjusted_p <- stats::p.adjust(res$adjusted_p_val, method = correction_method)
  res[, "adjusted_p_val" := adjusted_p]
  
  significant_indeces <- which(res$adjusted_p_val <= significant)
  if (length(significant_indeces) == 0) {
    warning("No significant terms were found")
    return()
  }
  
  if (contribution) {
    sig_cols <- columnSignificance(scores, gmt, background, cutoff,
                                   significant, correction_method, res$adjusted_p_val)
    res <- cbind(res, sig_cols[, -1])
  } else {
    sig_cols <- NULL
  }
  
  # if significant result were found and cytoscape file tag exists
  # proceed with writing files in the working directory
  if (length(significant_indeces) > 0 & !is.na(cytoscape_file_tag)) {
    prepareCytoscape(res[significant_indeces, c("term_id", "term_name", "adjusted_p_val")],
                     gmt[significant_indeces], 
                     cytoscape_file_tag,
                     sig_cols[significant_indeces,], color_palette, custom_colors, color_integrated_only)
  }
  
  res[significant_indeces]
}


#' Perform pathway enrichment analysis on an ordered list of genes
#'
#' @param genelist character vector of gene names, in decreasing order
#'   of significance
#' @param gmt GMT object
#' @param background character vector of gene names. List of all genes being used
#'   as a statistical background
#'
#' @return a data.table of terms with the following columns:
#'   \describe{
#'     \item{term_id}{The id of the term}
#'     \item{term_name}{The full name of the term}
#'     \item{adjusted_p_val}{The associated p-value adjusted for multiple testing}
#'     \item{term_size}{The number of genes annotated to the term}
#'     \item{overlap}{A character vector of the genes that overlap between the
#'        term and the query}
#'   }
#' @keywords internal
enrichmentAnalysis <- function(genelist, gmt, background) {
  dt <- data.table(term_id=names(gmt))
  
  for (i in 1:length(gmt)) {
    term <- gmt[[i]]
    tmp <- orderedHypergeometric(genelist, background, term$genes)
    overlap <- genelist[1:tmp$ind]
    overlap <- overlap[overlap %in% term$genes]
    if (length(overlap) == 0) overlap <- c()
    set(dt, i, 'term_name', term$name)
    set(dt, i, 'adjusted_p_val', tmp$p_val)
    set(dt, i, 'term_size', length(term$genes))
    set(dt, i, 'overlap', list(list(overlap)))
  }
  dt
}

#' Determine which terms are found to be significant using each column
#' individually. 
#'
#' @inheritParams ActivePathways
#' @param pvals p-value for the pathways calculated by ActivePathways
#'
#' @return a data.table with columns 'term_id' and a column for each column
#' in \code{scores}, indicating whether each term (pathway) was found to be
#' significant or not when considering only that column. For each term, 
#' either report the list of related genes if that term was significant, or NA if not. 

columnSignificance <- function(scores, gmt, background, cutoff, significant, correction_method, pvals) {
  dt <- data.table(term_id=names(gmt), evidence=NA)
  for (col in colnames(scores)) {
    col_scores <- scores[, col, drop=TRUE]
    col_scores <- col_scores[col_scores <= cutoff]
    col_scores <- names(col_scores)[order(col_scores)]
    
    res <- enrichmentAnalysis(col_scores, gmt, background)
    set(res, i = NULL, "adjusted_p_val", stats::p.adjust(res$adjusted_p_val, correction_method))
    set(res, i = which(res$adjusted_p_val > significant), "overlap", list(list(NA)))
    set(dt, i=NULL, col, res$overlap)
  }
  
  ev_names = colnames(dt[,-1:-2])
  set_evidence <- function(x) {
    ev <- ev_names[!is.na(dt[x, -1:-2])]
    if(length(ev) == 0) {
      if (pvals[x] <= significant) {
        ev <- 'combined'
      } else {
        ev <- 'none'
      }
    }
    ev
  }
  evidence <- lapply(1:nrow(dt), set_evidence)
  
  set(dt, i=NULL, "evidence", evidence)
  colnames(dt)[-1:-2] = paste0("Genes_", colnames(dt)[-1:-2])
  
  dt
}

#' Export the results from ActivePathways as a comma-separated values (CSV) file. 
#'
#' @param res the data.table object with ActivePathways results.
#' @param file_name location and name of the CSV file to write to.
#' @export
#'
#' @examples
#'     fname_scores <- system.file("extdata", "Adenocarcinoma_scores_subset.tsv", 
#'          package = "ActivePathways")
#'     fname_GMT = system.file("extdata", "hsapiens_REAC_subset.gmt",
#'          package = "ActivePathways")
#'
#'     dat <- as.matrix(read.table(fname_scores, header = TRUE, row.names = 'Gene'))
#'     dat[is.na(dat)] <- 1
#'
#'     res <- ActivePathways(dat, fname_GMT)
#'\donttest{
#'     export_as_CSV(res, "results_ActivePathways.csv")
#'}
export_as_CSV = function (res, file_name) {
	data.table::fwrite(res, file_name)	
} 

