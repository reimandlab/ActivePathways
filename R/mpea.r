#' mpea
#'
#' @param scores A numerical matrix of p-values where each row is a gene and
#'   each column is a test. Rownames should be the genes and colnames the names
#'   of the tests. All values must be 0<=p<=1 with missing values removed or
#'   converted to 1
#' @param gmt A GMT object to be used for enrichment analysis. If a filename, a
#'   GMT object will be read from the file
#' @param merge.method Method to merge p-values. See section Merging p Values
#' @param correction.method Method to correct p-values. See
#'   \code{\link[stats]{p.adjust}} for details
#' @param background A character vector of gene names to be used as a
#'   statistical background. By default, the background is all genes that appear
#'   in \code{gmt}
#' @param cutoff A maximum p-value for a gene to be used for enrichment analysis.
#'   Any genes with \code{p.val > significant} will be discarded before testing
#' @param significant A number in [0,1] denoting the maximum p-value for a
#'   pathway to be considered significantly enriched.
#' @param return.all Whether to return results for all terms or only significant
#'   terms
#' @param contribution Evaluate contribution of individual columns to a term's
#'   significance. See section Column Contribution
#' @param cytoscape.filenames a vector of 2 or 3 filenames denoting where
#'   information for Cytoscape will be written to (see section Cytoscape). If
#'   NULL, do not write any files.
#'
#' @return A data.table of terms containing the following columns:
#'   \describe{
#'     \item{term.id}{The id of the term}
#'     \item{term.name}{The full name of the term}
#'     \item{p.val}{The associated p-value}
#'     \item{term.size}{The number of genes annotated to the term}
#'     \item{overlap}{A character vector of the genes that overlap between the
#'     term and the query}
#'     \item{...}{If \code{contribution == TRUE}, a column for each test in
#'       \code{scores} denoting the contribution of that column to each
#'       pathway's p-value. See secion Column Contribution.}
#'   }
#'   If \code{return.all == FALSE} then only terms with
#'     \code{p.val <= significant} will be returned, otherwise all terms will be
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
#' @section Column Contribution: If \code{contribution == TRUE}, mpea will find
#'   the contribution of each column in \code{scores} to the p-value for each
#'   pathway. The contribution is reported as the log-fold-change in p-values
#'   when mpea is run with the column excluded from the data
#'   (\code{-log10(p_with_column / p_without_column)}). A contribution of 0
#'   means the p-value is the same with or without that column, while a large
#'   positive number means the p-value is significantly smaller when that column
#'   is included in the data.
#'
#' @section Cytoscape:
#'   If \code{cytoscape.filenames} is supplied, \code{mpea} will write three
#'   files that can be used to build a network using Cytoscape and the
#'   EnrichmentMap and enhancedGraphics apps. The three fies written are:
#'   \describe{
#'     \item{cytoscape.filenames[1]}{A list of significant terms and the
#'     associated p-value. Only terms with \code{p.val <= significant} are
#'     written to this file}
#'     \item{cytoscape.filenames[2]}{A matrix indicating whether the significant
#'     pathways are found to be significant when considering only one column from
#'     \code{scores}. A 1 indicates that that term is significant using only that
#'     column to test for enrichment analysis}
#'     \item{cytoscape.filenames[3]}{A Shortened version of the supplied gmt
#'     file, containing only the terms in \code{cytoscape.filenames[1]}}
#'   }
#'   If \code{contribution == FALSE} the matrix of column significance will not
#'   be written. Only 2 file names need to be supplied, and if three are given
#'   the second will be ignored
#'
#'   How to use: Create an enrichment map in Cytoscape with the file of terms
#'   (cytoscape.filenames[1]) and the shortened gmt file
#'   (cytoscape.filenames[3]). Upload (File > import > table > file) the
#'   subgroups file (cytoscape.filenames[2]) as a table. Under the 'style'
#'   panel, set image/Chart1 to use the column `instruct` and the passthrough
#'   mapping type.
#'
#' @examples
#' \dontrun{
#'     dat <- as.matrix(read.table('path/to/data.txt', header=TRUE, row.names='Gene'))
#'     dat[is.na(dat)] <- 1
#'     gmt <- read.GMT('path/to/gmt.gmt')
#'     mpea(dat, gmt, return.all=TRUE,
#'          cytoscape.filenames=c('terms.txt', 'groups.txt', 'abridged.gmt'))
#' }
#'
#' @export
#
# TODO: enter citations for article on merging p-values
# http://www.jstor.org/stable/2529826
# TODO: enter citations for Cytoscape, enrichmentMap, and enhancedGraphics
mpea <- function(scores, gmt, cutoff=0.1, significant=0.05, return.all=FALSE,
                 merge.method=c("Fisher", "Brown", "logitp", "meanp", "sump",
                                "sumz", "sumlog", "votep", "wilkinsonp"),
                 correction.method=c("holm", "fdr", "hochberg", "hommel",
                                     "bonferroni", "BH", "BY", "none"),
                 background=makeBackground(gmt), contribution=TRUE,
                 cytoscape.filenames=NULL) {


    merge.method <- match.arg(merge.method)
    correction.method <- match.arg(correction.method)

    ##### Validation #####

    if (!(is.matrix(scores) && is.numeric(scores))) stop("scores must be a numeric matrix")
    if (!is.numeric(cutoff) || cutoff < 0 || cutoff > 1) {
        stop("cutoff must be a value in [0,1]")
    }
    if (!is.numeric(significant) || significant < 0 || significant > 1) {
        stop("significant must be a value in [0,1]")
    }
    if (!is.character(background)) stop("background must be a character vector")
    if (!is.GMT(gmt)) gmt <- read.GMT(gmt)

    if (ncol(scores) == 1 && contribution) {
        contribution <- FALSE
        message("scores contains only one column. Column contributions will not be calculated")
    }
    if (!is.null(cytoscape.filenames)){
        if (contribution == TRUE && length(cytoscape.filenames) != 3) {
            stop("Must supply 3 file names to cytoscape.filenames")
        }
        if (!contribution){
            if (!length(cytoscape.filenames) %in% c(2,3)) {
                stop("Must supply 2 file names to cytoscape.filenames")
            }
            if (length(cytoscape.filenames) == 3) {
                message(paste("Column contributions will not be evaluated so the",
                              "contribution matrix is not being written.",
                              "cytoscape.filenames[2] will be ignored"))
                cytoscape.filenames <- cytoscape.filenames[-2]
            }
        }
    }

    ##### filtering and sorting #####

    # Remove any genes not found in the GMT
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
    res[, p.val := p.adjust(p.val, method=correction.method)]

    significant.indeces <- which(res$p.val <= significant)
    if (length(significant.indeces) == 0) {
        warning("No significant terms were found")
        if (!is.null(cytoscape.filenames)) warning("Cytoscape files were not written")
    }

    # If return.all==FALSE replace res and gmt with only the terms that will be
    # returned with the results
    if (!return.all) {
        res <- res[significant.indeces]
        gmt <- gmt[significant.indeces]
        significant.indeces <- 1:length(significant.indeces)
    }

    if (contribution) {
        col.contribution <- columnContribution(scores, gmt, background, cutoff,
                                               correction.method, merge.method, res$p.val)
        res <- merge(res, col.contribution, by='term.id', all=TRUE, sort=FALSE)
    }

    if (!is.null(cytoscape.filenames) && length(significant.indeces) > 0) {
        if (contribution) {
            sig.cols <- columnSignificance(scores, gmt[significant.indeces],
                                           background, cutoff, significant, correction.method)
        } else {
            sig.cols <- NULL
        }
        prepareCytoscape(res[significant.indeces, .(term.id, term.name, p.val)],
                         gmt[significant.indeces], cytoscape.filenames, sig.cols)
    }

    res
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
#'     \item{p.val}{The associated p-value}
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
        set(dt, i, 'p.val', tmp$p.val)
        set(dt, i, 'term.size', length(term$genes))
        set(dt, i, 'overlap', list(list(overlap)))
    }
    dt
}

#' Get contribution of each column to each term
#'
#' @param mpea.pvals p-values determined by mpea
#'
#' @inheritParams mpea
#'
#' @return a data.table of terms. The first column contains the term id. The
#'  rest of the columns give the log-fold-change when the column is excluded
#'
#' @examples
#' \dontrun{
#'   dat <- as.matrix(read.table('path/to/data.txt', header=TRUE, row.names='Gene'))
#'   dat[is.na(dat)] <- 1
#'   gmt <- read.GMT('path/to/gmt.gmt')
#'   res <- mpea(dat, gmt, contribution=FALSE)
#'   columnContribution(dat, gmt, makeBackground(gmt), 0.1, 'fdr', 'Fisher', res$p.val)
#' }
columnContribution <- function(scores, gmt, background, cutoff,
                               correction.method, merge.method, mpea.pvals) {
    dt <- data.table(term.id=names(gmt))

    for (col in colnames(scores)) {
        # Remove column and merge scores
        scores2 <- scores[, -which(col == colnames(scores)), drop=FALSE]
        merged.scores <- merge_p_values(scores2, merge.method)
        merged.scores <- merged.scores[merged.scores <= cutoff]
        merged.scores <- names(merged.scores)[order(merged.scores)]

        # Get p-valeus and report log-fold-change
        p.vals <- sapply(gmt, function(x)
            orderedHypergeometric(merged.scores, background, x$genes)$p.val)
        p.vals <- p.adjust(p.vals, method=correction.method)
        dt[, (col) := -log10(mpea.pvals / p.vals)]
     }
    dt
}

#' Determine which pathways are found to be significant using each column
#' individually
#'
#' @inheritParams mpea
#'
#' @return a data.table with columns 'term.id' and a column for each column
#' in \code{scores}, indicating whether each pathway was found to be
#' significant(TRUE) or not(FALSE) when considering only that column

columnSignificance <- function(scores, gmt, background, cutoff, significant, correction.method) {
    dt <- data.table(term.id=names(gmt))

    for (col in colnames(scores)) {
        col.scores <- scores[, col, drop=TRUE]
        col.scores <- col.scores[col.scores <= cutoff]
        col.scores <- names(col.scores)[order(col.scores)]

        p.vals <- sapply(gmt, function(x)
            orderedHypergeometric(col.scores, background, x$genes)$p.val)
        p.vals <- p.adjust(p.vals, correction.method)
        dt[, (col) := as.numeric(p.vals <= significant)]
    }
    dt
}
