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
#'       that column to test for enrichment analysis}. Left out if
#'       contribution==FALSE and therefore col.significance==NULL
#'     \item{cytoscape.filenames[3]}{A Shortened version of the supplied gmt
#'       file, containing only the terms in \code{cytoscape.filenames[1]}}
#'   }
#'
#' @param terms a data.table with columns 'term.id', 'term.name', 'p.val'
#' @param gmt an abridged gmt object containing only the pathways that were
#' found to be significant
#' @param filenames vector of 2 or 3 filesnames denoting where to write fles to
#' @param col.significance a data.table with a column 'term.id' and a column
#' for each test indicating whether a pathways is signficiant (TRUE) or not
#' (FALSE) when considering only that column. If contribution==TRUE, use
#' col.significance=NULL and this will be skipped
#'
#' @return None

prepareCytoscape <- function(terms, gmt, filenames, col.significance) {
    if (!is.null(col.significance)) {
        tests <- colnames(col.significance)[-1]
        rows <- 1:nrow(col.significance)

        ## Create groups file for enhancedgraphics. Keep both methods for now
        if (FALSE) {
            # Use outer ring
            col.significance[, contributes := paste(.SDcol, collapse='|'), by=rows, .SDcols=-1]
            colors <- paste(rainbow(length(tests), start=0.2), collapse=',')
            instruct.str <- paste('circoschart: firstarc=.7 arcwidth=.3 attributelist="contributes"',
                                  ' colorlist="', colors, '"',
                                  ' showlabels=FALSE', sep="")
            col.significance[, instruct := instruct.str]
        } else {
            # Use pichart
            col.significance[, none := as.numeric(all(!.SD)), by=rows, .SDcols=-1]
            colors <- paste(rainbow(length(tests)), collapse=',')
            instruct.str <- paste('piechart:',
                                  ' attributelist="', paste(tests, collapse=','), ',none"',
                                  ' colorlist="', colors, ',#CCCCCC"',
                                  ' showlabels=FALSE', sep='')
            col.significance[, instruct := instruct.str]
        }

        write.table(terms, file=filenames[1], row.names=FALSE, sep="\t", quote=FALSE)
        write.table(col.significance, file=filenames[2], row.names=FALSE, sep="\t", quote=FALSE)
        write.GMT(gmt, filenames[3])
    } else {
        write.table(terms, file=filenames[1], row.names=FALSE, sep="\t", quote=FALSE)
        write.GMT(gmt, filenames[2])
    }
}



