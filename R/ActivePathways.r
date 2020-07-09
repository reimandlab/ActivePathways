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
#' @param geneset.filter A numeric vector of length two giving the lower and 
#'   upper limits for the size of the annotated geneset to pathways in gmt.
#'   Pathways with a geneset shorter than \code{geneset.filter[1]} or longer
#'   than \code{geneset.filter[2]} will be removed. Set either value to NA to
#'   to not enforce a minimum or maximum value, or set \code{geneset.filter} to 
#'   \code{NULL} to skip filtering.
#' @param cutoff A maximum merged p-value for a gene to be used for analysis.
#'   Any genes with merged, unadjusted \code{p > significant} will be discarded 
#'   before testing.
#' @param significant Significance cutoff for selecting enriched pathways. Pathways with
#'   \code{adjusted.p.val < significant} will be selected as results.
#' @param merge.method Statistical method to merge p-values. See section on Merging P-Values
#' @param correction.method Statistical method to correct p-values. See
#'   \code{\link[stats]{p.adjust}} for details.
#' @param cytoscape.file.tag The directory and/or file prefix to which the output files
#'   for generating enrichment maps should be written. If NA, files will not be written. 
#'
#' @return A data.table of terms (enriched pathways) containing the following columns:
#'   \describe{
#'     \item{term.id}{The database ID of the term}
#'     \item{term.name}{The full name of the term}
#'     \item{adjusted.p.val}{The associated p-value, adjusted for multiple testing}
#'     \item{term.size}{The number of genes annotated to the term}
#'     \item{overlap}{A character vector of the genes enriched in the term}
#'     \item{evidence}{Columns of \code{scores} (i.e., omics datasets) that contributed 
#'          individually to the enrichment of the term. Each input column is evaluated 
#'          separately for enrichments and added to the evidence if the term is found.}
#'   }
#'
#' @section Merging P-values:
#' To obtain a single p-value for each gene across the multiple omics datasets considered, 
#' the p-values in \code{scores} #' are merged row-wise using a data fusion approach of p-value merging. 
#' The two available methods are:
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
#' }
#'
#' @section Cytoscape:
#'   To visualize and interpret enriched pathways, ActivePathways provides an option
#'   to further analyse results as enrichment maps in the Cytoscape software. 
#'   If \code{!is.na(cytoscape.file.tag)}, four files will be written that can be used 
#'   to build enrichment maps. This requires the EnrichmentMap and enhancedGraphics apps.
#'
#' The four files written are:
#'   \describe{
#'     \item{pathways.txt}{A list of significant terms and the
#'     associated p-value. Only terms with \code{adjusted.p.val <= significant} are
#'     written to this file.}
#'     \item{subgroups.txt}{A matrix indicating whether the significant terms (pathways)
#'     were also found to be significant when considering only one column from
#'     \code{scores}. A one indicates that that term was found to be significant 
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
                            geneset.filter = c(5, 1000), cutoff = 0.1, significant = 0.05,
                            merge.method = c("Brown", "Fisher"),
                            correction.method = c("holm", "fdr", "hochberg", "hommel",
                                                  "bonferroni", "BH", "BY", "none"),
                            cytoscape.file.tag = NA) {
  
  merge.method <- match.arg(merge.method)
  correction.method <- match.arg(correction.method)
  
  ##### Validation #####
  # scores
  if (!(is.matrix(scores) && is.numeric(scores))) stop("scores must be a numeric matrix")
  if (any(is.na(scores))) stop("scores may not contain missing values")
  if (any(scores < 0) || any(scores > 1)) stop("All values in scores must be in [0,1]")
  if (any(duplicated(rownames(scores)))) stop("Scores matrix contains duplicated genes - rownames must be unique.")
  
  # cutoff and significant
  stopifnot(length(cutoff) == 1)
  stopifnot(is.numeric(cutoff))
  if (cutoff < 0 || cutoff > 1) stop("cutoff must be a value in [0,1]")
  stopifnot(length(significant) == 1)
  stopifnot(is.numeric(significant))
  if (significant < 0 || significant > 1) stop("significant must be a value in [0,1]")
  
  # gmt
  if (!is.GMT(gmt)) gmt <- read.GMT(gmt)
  if (length(gmt) == 0) stop("No pathways in gmt made the geneset.filter")
  if (!(is.character(background) && is.vector(background))) {
    stop("background must be a character vector")
  } 
  
  # geneset.filter
  if (!is.null(geneset.filter)) {
    if (!(is.numeric(geneset.filter) && is.vector(geneset.filter))) {
      stop("geneset.filter must be a numeric vector")
    }
    if (length(geneset.filter) != 2) stop("geneset.filter must be length 2")
    if (!is.numeric(geneset.filter)) stop("geneset.filter must be numeric")
    if (any(geneset.filter < 0, na.rm=TRUE)) stop("geneset.filter limits must be positive")
  }
  
  contribution <- TRUE
  # contribution
  if (ncol(scores) == 1) {
    contribution <- FALSE
    message("Scores matrix contains only one column. Column contributions will not be calculated.")
  }
  
  ##### filtering and sorting ####
  
  # Filter the GMT
  if(!is.null(geneset.filter)) {
    orig.length <- length(gmt)
    if (!is.na(geneset.filter[1])) {
      gmt <- Filter(function(x) length(x$genes) >= geneset.filter[1], gmt)
    }
    if (!is.na(geneset.filter[2])) {
      gmt <- Filter(function(x) length(x$genes) <= geneset.filter[2], gmt)
    }
    if (length(gmt) == 0) stop("No pathways in gmt made the geneset.filter")
    if (length(gmt) < orig.length) {
      message(paste(orig.length - length(gmt), "terms were removed from gmt", 
                    "because they did not make the geneset.filter"))
    }
  }
  
  # Remove any genes not found in the background
  orig.length <- nrow(scores)
  scores <- scores[rownames(scores) %in% background, , drop=FALSE]
  if (nrow(scores) == 0) {
    stop("scores does not contain any genes in the background")
  }
  if (nrow(scores) < orig.length) {
    message(paste(orig.length - nrow(scores), "rows were removed from scores",
                  "because they are not found in the background"))
  }
  
  # merge p-values to get a single score for each gene and remove any genes
  # that don't make the cutoff
  merged.scores <- merge_p_values(scores, merge.method)
  merged.scores <- merged.scores[merged.scores <= cutoff]
  
  if (length(merged.scores) == 0) stop("No genes made the cutoff")
  
  # Sort genes by p-value
  ordered.scores <- names(merged.scores)[order(merged.scores)]
  
  ##### enrichmentAnalysis and column contribution #####
  
  res <- enrichmentAnalysis(ordered.scores, gmt, background)
  adjusted_p <- stats::p.adjust(res$adjusted.p.val, method = correction.method)
  res[, "adjusted.p.val" := adjusted_p]
  
  significant.indeces <- which(res$adjusted.p.val <= significant)
  if (length(significant.indeces) == 0) {
    warning("No significant terms were found.")
    return()
  }
  
  if (contribution) {
    sig.cols <- columnSignificance(scores, gmt, background, cutoff,
                                   significant, correction.method, res$adjusted.p.val)
    res <- cbind(res, sig.cols[, -1])
  } else {
    sig.cols <- NULL
  }
  
  # if significant result were found and cytoscape file tag exists
  # proceed with writing files in the working directory
  if (length(significant.indeces) > 0 & !is.na(cytoscape.file.tag)) {
    prepareCytoscape(res[significant.indeces, c("term.id", "term.name", "adjusted.p.val")],
                     gmt[significant.indeces], 
                     cytoscape.file.tag,
                     sig.cols[significant.indeces,])
  }
  
  res[significant.indeces]
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
#'     \item{term.id}{The id of the term}
#'     \item{term.name}{The full name of the term}
#'     \item{adjusted.p.val}{The associated p-value adjusted for multiple testing}
#'     \item{term.size}{The number of genes annotated to the term}
#'     \item{overlap}{A character vector of the genes that overlap between the
#'        term and the query}
#'   }
#' @keywords internal
enrichmentAnalysis <- function(genelist, gmt, background) {
  dt <- data.table(term.id=names(gmt))
  
  for (i in 1:length(gmt)) {
    term <- gmt[[i]]
    tmp <- orderedHypergeometric(genelist, background, term$genes)
    overlap <- genelist[1:tmp$ind]
    overlap <- overlap[overlap %in% term$genes]
    if (length(overlap) == 0) overlap <- NA
    set(dt, i, 'term.name', term$name)
    set(dt, i, 'adjusted.p.val', tmp$p.val)
    set(dt, i, 'term.size', length(term$genes))
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
#' @return a data.table with columns 'term.id' and a column for each column
#' in \code{scores}, indicating whether each term (pathway) was found to be
#' significant or not when considering only that column. For each term, 
#' either report the list of related genes if that term was significant, or NA if not. 

columnSignificance <- function(scores, gmt, background, cutoff, significant, correction.method, pvals) {
  dt <- data.table(term.id=names(gmt), evidence=NA)
  for (col in colnames(scores)) {
    col.scores <- scores[, col, drop=TRUE]
    col.scores <- col.scores[col.scores <= cutoff]
    col.scores <- names(col.scores)[order(col.scores)]
    
    res <- enrichmentAnalysis(col.scores, gmt, background)
    set(res, i = NULL, "adjusted.p.val", stats::p.adjust(res$adjusted.p.val, correction.method))
    set(res, i = which(res$adjusted.p.val > significant), "overlap", list(list(NA)))
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

