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
#' @param filenames vector of 2 or 3 filesnames denoting where to write fles to
#'
#' @inheritParams mpea
#'
#' @return None
#'
#' @examples
#' \dontrun{
#'     res <- mpea(dat, gmt)
#'     prepareCytoscape(res, gmt, paste('cytoscape', c('terms.txt', 'groups.txt', 'abridged.gmt'), sep='/'), TRUE)
#' }

prepareCytoscape <- function(terms, gmt, filenames, contribution) {
    termlist <- terms[, .(term.id, term.name, p.val)]

    abridged.gmt <- gmt[unlist(terms$term.id)]

    if (contribution) {
        tests <- colnames(terms)[-1:-5]

        ## Create groups file for enhancedgraphics. Keep both methods for now
        if (FALSE) {
            # Use outer ring
            subgroups <- terms[, .(term.id)]
            contributes <- apply(terms, 1, function(x) paste(x[tests], collapse='|'))
            colors <- paste(rainbow(length(tests), start=0.2), collapse=',')
            instruct.str <- paste('circoschart: firstarc=.7 arcwidth=.3 attributelist="contributes"',
                                  ' colorlist="', colors, '"',
                                  ' showlabels=FALSE', sep="")
            subgroups[, c('contributes', 'instruct') := list(contributes, instruct.str)]
        } else {
            # Use pichart
            subgroups <- terms[, c('term.id', tests), with=FALSE]
            subgroups[, none := as.numeric(all(!.SD)), by=1:nrow(subgroups), .SDcols=tests]
            colors <- paste(rainbow(length(tests)), collapse=',')
            instruct.str <- paste('piechart:',
                                  ' attributelist="', paste(tests, collapse=','), ',none"',
                                  ' colorlist="', colors, ',#CCCCCC"',
                                  ' showlabels=FALSE', sep='')
            subgroups[, 'instruct' := instruct.str]
        }

        write.table(termlist, file=filenames[1], row.names=FALSE, sep="\t", quote=FALSE)
        write.table(subgroups, file=filenames[2], row.names=FALSE, sep="\t", quote=FALSE)
        write.GMT(abridged.gmt, filenames[3])
    } else {
        write.table(termlist, file=filenames[1], row.names=FALSE, sep="\t", quote=FALSE)
        write.GMT(abridged.gmt, filenames[2])
    }
}



