#' Prepare files for building an enrichment map network visualization in Cytoscape
#'
#' This function writes four text files that are used to build an network using
#' Cytoscape and the EnrichmentMap app. The files are prefixed with \code{cytoscape.file.tag}. 
#'   The four files written are:
#'   \describe{
#'     \item{pathways.txt}{A list of significant terms and the
#'     associated p-value. Only terms with \code{adjusted.p.val <= significant} are
#'     written to this file}
#'     \item{subgroups.txt}{A matrix indicating whether the significant
#'     pathways are found to be significant when considering only one column (i.e., type of omics evidence) from
#'     \code{scores}. A 1 indicates that that term is significant using only that
#'     column to test for enrichment analysis}
#'     \item{pathways.gmt}{A shortened version of the supplied GMT
#'     file, containing only the terms in pathways.txt.}
#'     \item{legend.pdf}{A legend with colours matching contributions
#'     from columns in \code{scores}}
#'   }
#'
#' @param terms A data.table object with the columns 'term.id', 'term.name', 'adjusted.p.val'. 
#' @param gmt An abridged GMT object containing only the pathways that were
#' found to be significant in the ActivePathways analysis.
#' @param cytoscape.file.tag The user-defined file prefix and/or directory defining the location of the files.
#' @param col.significance A data.table object with a column 'term.id' and a column
#' for each type of omics evidence indicating whether a term was also found to be signficiant or not
#' when considering only the genes and p-values in the corresponding column of the \code{scores} matrix. If term was not found, NA's are shown in columns, otherwise the relevant lists of genes are shown.
#' @import ggplot2
#'
#' @return None

prepareCytoscape <- function(terms, 
                             gmt, 
                             cytoscape.file.tag, 
                             col.significance) {
  if (!is.null(col.significance)) {
    tests <- unique(unlist(col.significance$evidence))
    rows <- 1:nrow(col.significance)
    
    evidence.columns = do.call(rbind, lapply(col.significance$evidence,
                                             function(x) 0+(tests %in% x)))
    colnames(evidence.columns) = tests
    
    col.significance = cbind(col.significance[,"term.id"], evidence.columns)
    
    # Use pichart
    col.colors <- grDevices::rainbow(length(tests))
    instruct.str <- paste('piechart:',
                          ' attributelist="', 
                          paste(tests, collapse=','),
                          '" colorlist="', 
                          paste(col.colors, collapse=','), 
                          '" showlabels=FALSE', sep='')
    col.significance[, "instruct" := instruct.str]
    
    # Writing the Files
    utils::write.table(terms, 
                file=paste0(cytoscape.file.tag, "pathways.txt"), 
                row.names=FALSE, 
                sep="\t", 
                quote=FALSE)
    utils::write.table(col.significance, 
                file=paste0(cytoscape.file.tag, "subgroups.txt"), 
                row.names=FALSE, 
                sep="\t", 
                quote=FALSE)
    write.GMT(gmt, 
              paste0(cytoscape.file.tag, "pathways.gmt"))
    
    # Making a Legend
      dummy_plot = ggplot(data.frame("tests" = factor(tests, levels = tests),
                                     "value" = 1), aes(tests, fill = tests)) +
        geom_bar() +
        scale_fill_manual(name = "Contribution", values=col.colors)

      grDevices::pdf(file = NULL) # Suppressing Blank Display Device from ggplot_gtable
      dummy_table = ggplot_gtable(ggplot_build(dummy_plot))
      grDevices::dev.off()

      legend = dummy_table$grobs[[which(sapply(dummy_table$grobs, function(x) x$name) == "guide-box")]]
      
      # Estimating height & width
      legend_height = ifelse(length(tests) > 20, 
                             5.5, 
                             length(tests)*0.25+1)
      legend_width = ifelse(length(tests) > 20, 
                            ceiling(length(tests)/20)*(max(nchar(tests))*0.05+1), 
                            max(nchar(tests))*0.05+1)
      ggsave(legend,
             device = "pdf",
             filename = paste0(cytoscape.file.tag, "legend.pdf"), 
             height = legend_height, 
             width = legend_width, 
             scale = 1)
    
  } else {
    utils::write.table(terms, 
                file=paste0(cytoscape.file.tag, "pathways.txt"),
                row.names=FALSE, 
                sep="\t", 
                quote=FALSE)
    write.GMT(gmt, 
              paste0(cytoscape.file.tag, "pathways.gmt"))
  }
}



