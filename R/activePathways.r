#' activePathways
#'
#' @param scores A numerical matrix of p-values where each row is a gene and
#'   each column is a test. Rownames should be the genes and colnames the names
#'   of the tests. All values must be 0<=p<=1 with missing values removed or
#'   converted to 1
#' @param gmt A GMT object to be used for enrichment analysis. If a filename, a
#'   GMT object will be read from the file
#' @param background A character vector of gene names to be used as a
#'   statistical background. By default, the background is all genes that appear
#'   in \code{gmt}
#' @param geneset.filter A numeric vector of length two giving the lower and 
#'   upper limits for the size of the annotated geneset to pathways in gmt.
#'   Pathways with a geneset shorter than \code{geneset.filter[1]} or longer
#'   than \code{geneset.filter[2]} will be removed. Set either value to NA to
#'   to not enforce a minimum or maximum value, or set \code{geneset.filter} to 
#'   \code{NULL} to skip filtering
#' @param cutoff A maximum p-value for a gene to be used for enrichment analysis.
#'   Any genes with \code{adjusted.p.val > significant} will be discarded before testing
#' @param significant A number in [0,1] denoting the maximum p-value for a
#'   pathway to be considered significantly enriched.
#' @param merge.method Method to merge p-values. See section Merging p Values
#' @param correction.method Method to correct p-values. See
#'   \code{\link[stats]{p.adjust}} for details
#' @param return.all Whether to return results for all terms or only significant
#'   terms
#' @param cytoscape.file.tag a user-chosen tag to identify this instance of ActivePathways and
#'   used to name cytoscape output files. If NULL, do not write any files.
#' @param cytoscape.file.dir the directory to which the output files should be written, if unspecified,
#'   files will be written to the current directory. If directory does not exist, it will be automatically
#'   created.
#'
#' @return A data.table of terms containing the following columns:
#'   \describe{
#'     \item{term.id}{The id of the term}
#'     \item{term.name}{The full name of the term}
#'     \item{adjusted.p.val}{The associated p-value, adjusted for multiple testing}
#'     \item{term.size}{The number of genes annotated to the term}
#'     \item{overlap}{A character vector of the genes that overlap between the
#'     term and the query}
#'     \item{evidence}{Columns of \code{scores} that contributed individually to 
#'          enrichment of the pathway. Each column is evaluated separately for 
#'          enrichments and added to the evidence field if the pathway is found.}
#'   }
#'   If \code{return.all == FALSE} then only terms with
#'     \code{adjusted.p.val <= significant} will be returned, otherwise all terms will be
#'     returned.
#'
#' @section Merging p Values:
#' In order to obtain a single score for each gene, the p-values in \code{scores}
#' are merged row-wise. There are multiple methods available that can be used
#' to obtain this merged score. The main methods are:
#' \describe{
#'  \item{Fisher or sumlog}{Fisher's method assumes p-values are uniformly
#'  distributed and performs a chi-squared test on the statistic sum(-2 log(p)).
#'  This method is most appropriate when the columns in \code{scores} are
#'  independent.}
#'  \item{Brown}{Brown's method extends Fisher's method by accounting for the
#'  covariance in the columns of \code{scores}. It is more appropriate when the
#'  tests of significance used to create the columns in \code{scores} are not
#'  necessarily independent.}
#' }
#' Other methods are also available. See \code{\link[metap]{metap-package}}
#' for more details
#'
#' @section Cytoscape:
#'   If \code{cytoscape.file.tag} is supplied, activePathways will write four
#'   files that can be used to build a network using Cytoscape and the
#'   EnrichmentMap and enhancedGraphics apps. The four files written are:
#'   \describe{
#'     \item{\code{cytoscape.file.tag}_pathways.txt}{A list of significant terms and the
#'     associated p-value. Only terms with \code{adjusted.p.val <= significant} are
#'     written to this file}
#'     \item{\code{cytoscape.file.tag}_subgroups.txt}{A matrix indicating whether the significant
#'     pathways are found to be significant when considering only one column from
#'     \code{scores}. A 1 indicates that that term is significant using only that
#'     column to test for enrichment analysis}
#'     \item{\code{cytoscape.file.tag}.gmt}{A Shortened version of the supplied gmt
#'     file, containing only the terms in \code{cytoscape.file.tag}_pathways.txt}
#'     \item{\code{cytoscape.file.tag}_legend.pdf}{A legend with colours matching contributions
#'     from columns in \code{scores}}
#'   }
#'
#'   How to use: Create an enrichment map in Cytoscape with the file of terms
#'   (\code{cytoscape.file.tag}_pathways.txt) and the shortened gmt file
#'   (\code{cytoscape.file.tag}.gmt). Upload (File > import > table > file) the
#'   subgroups file (\code{cytoscape.file.tag}_subgroups.txt) as a table. Under the 'style'
#'   panel, set image/Chart1 to use the column `instruct` and the passthrough
#'   mapping type. Use \code{cytoscape.file.tag}_legend.pdf as a reference in final figure.
#'
#' @examples
#' \dontrun{
#'     dat <- as.matrix(read.table('path/to/data.txt', header=TRUE, row.names='Gene'))
#'     dat[is.na(dat)] <- 1
#'     activePathways(dat, 'path/to/gmt.gmt', return.all=TRUE,
#'          cytoscape.file.tag="dat",
#'          cytoscape.file.dir="results")
#' }
#'
#' @import data.table
#' @import metap
#'
#' @export

# TODO: enter citations for article on merging p-values
# http://www.jstor.org/stable/2529826
# TODO: enter citations for Cytoscape, enrichmentMap, and enhancedGraphics
activePathways <-  function(scores, gmt, background = makeBackground(gmt),
                            geneset.filter = c(5, 1000), cutoff = 0.1, significant = 0.05,
                            merge.method = c("Brown", "Fisher", "logitp", "meanp", 
                                             "sump", "sumz", "sumlog"),
                            correction.method = c("holm", "fdr", "hochberg", "hommel",
                                                  "bonferroni", "BH", "BY", "none"),
                            return.all=FALSE, 
                            cytoscape.file.tag = NULL,
                            cytoscape.file.dir = "") {
  
  merge.method <- match.arg(merge.method)
  correction.method <- match.arg(correction.method)
  
  ##### Validation #####
  # scores
  if (!(is.matrix(scores) && is.numeric(scores))) stop("scores must be a numeric matrix")
  if (any(is.na(scores))) stop("scores may not contain missing values")
  if (any(scores < 0) || any(scores > 1)) stop("All values in scores must be in [0,1]")
  if (any(duplicated(rownames(scores)))) stop("scores contains duplicated genes. rownames must be unique")
  
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
    message("scores contains only one column. Column contributions will not be calculated")
  }
  
  # cytoscape.file.dir
  if(!is.character(cytoscape.file.dir) | length(cytoscape.file.dir) != 1){
    stop("cytoscape.file.dir must be a string")
  }
  if(cytoscape.file.dir != "" & !is.null(cytoscape.file.tag)){
    if (!dir.exists(cytoscape.file.dir)){
      dir.create(cytoscape.file.dir)
      message(paste0("Creating ", cytoscape.file.dir))
      }
  }
  if(!endsWith(cytoscape.file.dir, "[/]")){
    cytoscape.file.dir = paste0(cytoscape.file.dir, "/")
  }
  
  # cytoscape.file.tag
  if(!is.character(cytoscape.file.tag) | length(cytoscape.file.tag) != 1){
    stop("cytoscape.file.tag must be a string")
  }
  
  # cytoscape.filenames # No longer in use
  # if (!is.null(cytoscape.filenames)){
  #   if (contribution == TRUE && length(cytoscape.filenames) != 4) {
  #     stop("Must supply 4 file names to cytoscape.filenames")
  #   }
  #   if (!contribution){
  #     if (!length(cytoscape.filenames) %in% c(3,4)) {
  #       stop("Must supply 3 file names to cytoscape.filenames")
  #     }
  #     if (length(cytoscape.filenames) == 4) {
  #       message(paste("Column contributions will not be evaluated so the",
  #                     "contribution matrix is not being written.",
  #                     "cytoscape.filenames[2] will be ignored"))
  #       cytoscape.filenames <- cytoscape.filenames[-2]
  #     }
  #   }
  # }
  
  ##### filtering and sorting #####
  
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
  res[, adjusted.p.val := p.adjust(adjusted.p.val, method=correction.method)]
  
  significant.indeces <- which(res$adjusted.p.val <= significant)
  if (length(significant.indeces) == 0) {
    warning("No significant terms were found")
    if (!is.null(cytoscape.file.tag)) warning("Cytoscape files were not written")
  }
  
  if (contribution) {
    sig.cols <- columnSignificance(scores, gmt, background, cutoff,
                                   significant, correction.method, res$adjusted.p.val)
    res <- cbind(res, sig.cols[, -1])
  } else {
    sig.cols <- NULL
  }
  
  if (!is.null(cytoscape.file.tag) && length(significant.indeces) > 0) {
    prepareCytoscape(res[significant.indeces, .(term.id, term.name, adjusted.p.val)],
                     gmt[significant.indeces], 
                     cytoscape.file.tag,
                     cytoscape.file.dir,
                     sig.cols[significant.indeces,])
  }
  
  if(return.all) return(res)
  res[significant.indeces]
}


#' Perform Gene Set Enrichment Analysis on an ordered list of genes
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
#'
#' @examples
#' \dontrun{
#'     enrichmentAnalysis(c('HERC2', 'SMC5', 'XPC', 'WRN'), gmt, makeBackground(gmt))
#' }
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

#' Determine which pathways are found to be significant using each column
#' individually
#'
#' @inheritParams activePathways
#' @param pvals p-value for the pathways calculated by activePathways
#'
#' @return a data.table with columns 'term.id' and a column for each column
#' in \code{scores}, indicating whether each pathway was found to be
#' significant(TRUE) or not(FALSE) when considering only that column

columnSignificance <- function(scores, gmt, background, cutoff, significant, correction.method, pvals) {
  dt <- data.table(term.id=names(gmt), evidence=NA)
  for (col in colnames(scores)) {
    col.scores <- scores[, col, drop=TRUE]
    col.scores <- col.scores[col.scores <= cutoff]
    col.scores <- names(col.scores)[order(col.scores)]
    
    res <- enrichmentAnalysis(col.scores, gmt, background)
    set(res, i=NULL, "adjusted.p.val", p.adjust(res$adjusted.p.val, correction.method))
    set(res, i=which(res$adjusted.p.val>significant), "overlap", list(list(NA)))
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

