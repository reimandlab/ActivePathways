#' Prepare files for building an Enrichment Map in Cytoscape
#'
#' This function writes three files that are used to build an network using
#' Cytoscape and the EnrichmentMap app.
#'   The four files written are:
#'   \describe{
#'     \item{cytoscape.filenames[1]}{A list of significant terms and the
#'       associated adjusted p-value. Only terms with \code{adjusted.p.val <= significant} are
#'       written to this file}
#'     \item{cytoscape.filenames[2]}{A matrix of column contributions to each
#'       term. A \code{1} indicates that that term is significant using only
#'       that column to test for enrichment analysis}. Left out if
#'       contribution==FALSE and therefore col.significance==NULL
#'     \item{cytoscape.filenames[3]}{A Shortened version of the supplied gmt
#'       file, containing only the terms in \code{cytoscape.filenames[1]}}
#'     \item{cytoscape.filenames[4]}{A legend with colours matching contributions
#'     from columns in \code{scores}}
#'   }
#'
#' @param terms a data.table with columns 'term.id', 'term.name', 'adjusted.p.val'
#' @param gmt an abridged gmt object containing only the pathways that were
#' found to be significant
#' @param filenames vector of 3 or 4 filesnames denoting where to write files to
#' @param col.significance a data.table with a column 'term.id' and a column
#' for each test indicating whether a pathway is signficiant (TRUE) or not
#' (FALSE) when considering only that column. If contribution==TRUE, use
#' col.significance=NULL and this will be skipped
#'
#' @return None

prepareCytoscape <- function(terms, gmt, filenames, col.significance) {
  if (!is.null(col.significance)) {
    tests <- unique(unlist(col.significance$evidence))
    rows <- 1:nrow(col.significance)
    
    evidence.columns = do.call(rbind, lapply(col.significance$evidence,
                                             function(x) 0+(tests %in% x)))
    colnames(evidence.columns) = tests
    
    col.significance = cbind(col.significance[,"term.id"], evidence.columns)
    
    # Use pichart
    col.colors <- paste(rainbow(length(tests)), collapse=',')
    instruct.str <- paste('piechart:',
                          ' attributelist="', 
                          paste(tests, collapse=','),
                          '" colorlist="', 
                          paste(col.colors, collapse=','), 
                          '" showlabels=FALSE', sep='')
    col.significance[, instruct := instruct.str]
    
    # Writing the Files
    write.table(terms, file=filenames[1], row.names=FALSE, sep="\t", quote=FALSE)
    write.table(col.significance, file=filenames[2], row.names=FALSE, sep="\t", quote=FALSE)
    write.GMT(gmt, filenames[3])
    
    # Making a Legend
    require(ggplot2)
    dummy_plot = ggplot(data.frame(test, 1), aes(test, fill = test)) +
      geom_bar() +
      scale_fill_manual(name = "Contribution", values=col.colors)
    dummy_table = ggplot_gtable(ggplot_build(dummy_plot))
    legend = dummy_table$grobs[[which(sapply(tmp$grobs, function(x) x$name) == "guide-box")]]
    legend_height = ifelse(length(test) > 19, 5, length(test)*0.2 +1)
    legend_width = ifelse(length(test) > 20, ceiling(length(test)/20)*(max(nchar(default_test))*0.05+1), max(nchar(default_test))*0.05+1)
    ggsave(legend, filename = filenames[4], height = legend_height, width = legend_width, scale = 1)
    
  } else {
    write.table(terms, file=filenames[1], row.names=FALSE, sep="\t", quote=FALSE)
    write.GMT(gmt, filenames[2])
  }
}



