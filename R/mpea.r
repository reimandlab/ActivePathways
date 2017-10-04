#' Title
#'
#' Summary
#'
#' Description
#'
#' @param scores A numerical matrix of p-values where each row is a gene and
#'   each column is a test. Rownames should be the genes and colnames the names
#'   of the tests. All values must be 0<=p<=1 with missing values converted to 1
#' @param gmt A GMT object to be used for enrichment analysis. If a filename,
#'   a GMT object will be created from the file
#' @param merge.method Method to merge p-values. See section Merging p-Values
#' @param correction.method Method to correct p-values. See Adjusting p-Values
#' @param background A character vector of gene names to be used as a
#'   statistical background. By default, the background is all genes that
#'   appear in \code{gmt}
#' @param cutoff A maximum p-value for agene to be used for enrichment analysis.
#'   Any genes with \code{p.val > significant} will be discarded before testing
#' @param significant: A number in [0,1] denoting the maximum p-value for a
#'   pathway to be considered significantly enriched.
#' @param return.all Whether to return results for all terms or only significant
#'   terms
#' @param contribution Evaluate contribution of individual columns to a term's
#'   significance. See section Column Contribution
#' @param cytoscape.filenames: a vector of 2 or 3 filenames denoting where information
#'   for Cytoscape will be written to (see section Cytoscape). If NULL,
#'   do not write any files.
#'
#' @return A data.table of terms containing the following columns:
#'   \describe{
#'     \item{term.id}{The id of the term}
#'     \item{term.name}{The full name of the term}
#'     \item{p.val}{The associated p-value}
#'     \item{term.size}{The number of genes annotated to the term}
#'     \item{overlap}{A character vector of the genes that overlap between the
#'        term and the query}
#'     \item{...}{If \core{contribution == TRUE}, A column for each test in
#'       \code{scores} denoting whether the term is found to be significant if
#'       using only that test in analysis}
#'   }
#'   If \code{return.all == FALSE} then only terms with
#'   \code{p.val <= significant} will be returned.
#'
#' @section Cytoscape:
#'   If cytoscape.filenames is supplied, mpea will write three files that can
#'   be used to build an network using Cytoscape and the EnrichmentMap app.
#'   The three fies written are:
#'   \describe{
#'     \item{cytoscape.filenames[1]}{A list of significant terms and the
#'       associated p-value. Only terms with \code{p.val <= significant} are
#'       written to this file}
#'     \item{cytoscape.filenames[2]}{A matrix of column contributions to each
#'       term. A \code{1} indicates that that term is significant using only
#'       that column to test for enrichment analysis}
#'     \item{cytoscape.filenames[3]}{A Shortened version of the supplied gmt
#'       file, containing only the terms in \code{cytoscape.filenames[1]}}
#'   }
#'   If \code{contribution == FALSE} the matrix of column contributions will
#'     not be written. Only 2 file names need to be supplied, and if three are
#'     given the second will be ignored
#'
#' @section Merging p-Values:
#'
#' @section Adjusting p-Values:
#'
#' @section Column Contribution:
#'   If \code{contribution == TRUE}, mpea will
#'
#' @export
mpea <- function(scores, gmt,
                 merge.method=c("Fisher", "Brown", "logtip", "meanp", "sump", "sumz", "sumlog"),
                 correction.method=c("fdr"), background=makeBackground(gmt),
                 cutoff=0.1, significant=0.05, return.all=FALSE, contribution=TRUE,
                 cytoscape.filenames=c('termlist.txt', 'subgroups.txt', 'abridged.gmt')) {

    # Validation
    merge.method <- match.arg(merge.method)
    correction.method <- match.arg(correction.method)

    if (!is.numeric(scores)) stop("scores must be a numeric matrix")
    if (!is.numeric(cutoff) || cutoff < 0 || cutoff > 1) stop("cutoff must be a value in [0,1]")
    if (!is.numeric(significant) || significant < 0 || significant > 1) {
        stop("significant must be a value in [0,1]")
    }
    if (!is.character(background)) stop("background must be a character vector")
    if (!is.GMT(gmt)) gmt <- read.GMT(gmt)

    if (!is.null(cytoscape.filenames)){
        if (contribution == TRUE && length(cytoscape.filenames) != 3) {
            stop("Must supply 3 file names to cytoscape.filenames")
        }
        if (contribution == FALSE){
            if (!length(cytoscape.filenames) %in% c(2,3)) {
                stop("Must supply 2 file names to cytoscape.filenames")
            }
            if (length(cytoscape.filenames) == 3) {
                warning("Column contributions will not be evaluated so the contribution matrix is not being written. cytoscape.filenames[2] will be ignored")
                cytoscape.filenames <- cytoscape.filenames[-2]
            }
        }
    }

    # Remove any genes not found in the GMT
    scores <- scores[rownames(scores) %in% background, , drop=FALSE]

    # merge p-values to get a single score for each gene and remove any genes
    # that don't make the cutoff
    merged.scores <- merge_p_values(scores, merge.method)
    merged.scores <- merged.scores[merged.scores <= cutoff]
    if (length(merged.scores) == 0) stop("No genes made the cutoff")

    # Sort genes by p-value
    ordered.scores <- names(merged.scores)[order(merged.scores)]

    res <- gsea(ordered.scores, gmt, background)
    res[, p.val := p.adjust(p.val, method=correction.method)]

    if (contribution){
        column.contribution <- columnContribution(scores, gmt, background, significant, cutoff, correction.method)
        res <- cbind(res, column.contribution[, -c('term.id', 'term.name')])
    }

    significant.indeces <- which(res$p.val <= significant)
    if (length(significant.indeces) == 0) {
        if (return.all) {
            warning("No significant terms were found")
        } else {
            stop("No significant terms were found")
        }
    } else {
        if (!is.null(cytoscape.filenames)) {
            prepareCytoscape(res[significant.indeces], gmt, cytoscape.filenames, contribution)
        }
    }

    if (return.all) return(res)
    return(res[significant.indeces])
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
gsea <- function(genelist, gmt, background) {
    dt <- data.table(term.id=names(gmt))

    for (i in 1:length(gmt)) {
        term <- gmt[[i]]
        tmp <- orderedHypergeometric(genelist, background, term$genes)
        overlap <- genelist[1:tmp$ind][genelist[1:tmp$ind] %in% term$genes]
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
#' @inheritParams mpea
#'
#' @return a data.table of terms. The first two columns contain the term id and
#'   name. The rest of the columns denoting whether the term is found to be
#'   significant if using only that test in analysis
columnContribution <- function(scores, gmt, background, significant, cutoff, correction.method) {
    dt <- data.table(term.id=names(gmt), term.name=sapply(gmt, function(x) x$name))

    for (col in colnames(scores)) {
        col.scores <- scores[, col, drop=TRUE]
        col.scores <- col.scores[col.scores <= cutoff]
        col.scores <- names(col.scores[order(col.scores)])

        p.vals <- sapply(gmt, function(x) orderedHypergeometric(col.scores, background, x$genes)$p.val)
        p.vals <- p.adjust(p.vals, method=correction.method)
        p.vals <- as.numeric(p.vals <= significant)
        dt[, (col) := p.vals]
     }
    dt
}

#' Prepare files for building an Enrichment Map in Cytoscape
#'
#' This function writes three files that are used to build an network using
#' Cytoscape and the EnrichmentMap app.
#'   The three fies written are:
#'   \describe{
#'     \item{cytoscape.filenames[1]}{A list of significant terms and the
#'       associated p-value. Only terms with \code{p.val <= significant} are
#'       written to this file}
#'     \item{cytoscape.filenames[2]}{A matrix of column contributions to each
#'       term. A \code{1} indicates that that term is significant using only
#'       that column to test for enrichment analysis}
#'     \item{cytoscape.filenames[3]}{A Shortened version of the supplied gmt
#'       file, containing only the terms in \code{cytoscape.filenames[1]}}
#'   }
#'
#' @param terms a data.table of terms, the output of mpea()
#'
#' @inheritParams mpea
#'
#' @return None

prepareCytoscape <- function(terms, gmt, filenames, contribution) {
    termlist <- terms[, .(term.id, term.name, p.val)]

    abridged.gmt <- gmt[unlist(terms$term.id)]

    if (contribution) {
        tests <- colnames(terms)
        tests <- tests[6:length(tests)]
        subgroups <- terms[, c('term.id', tests), with=FALSE]

        tests <- colnames(subgroups[, -'term.id'])
        instruct.str <- paste('pichart:',
                              ' attributelist="', paste(tests,collapse=","), '"',
                              ' colorlist="', paste(rainbow(length(tests)), collapse=","), '"',
                              ' showlabels=FALSE', sep="")
        subgroups[, instruct := instruct.str]
        write.table(termlist, file=filenames[1], row.names=FALSE, sep="\t", quote=FALSE)
        write.table(subgroups, file=filenames[2], row.names=FALSE, sep="\t", quote=FALSE)
        write.GMT(abridged.gmt, filenames[3])
    } else {
        write.table(termlist, file=filenames[1], row.names=FALSE, sep="\t", quote=FALSE)
        write.GMT(abridged.gmt, filenames[2])
    }



}



