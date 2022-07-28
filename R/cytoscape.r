#' Prepare files for building an enrichment map network visualization in Cytoscape
#'
#' This function writes four text files that are used to build an network using
#' Cytoscape and the EnrichmentMap app. The files are prefixed with \code{cytoscape_file_tag}. 
#'   The four files written are:
#'   \describe{
#'     \item{pathways.txt}{A list of significant terms and the
#'     associated p-value. Only terms with \code{adjusted_p_val <= significant} are
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
#' @param terms A data.table object with the columns 'term_id', 'term_name', 'adjusted_p_val'. 
#' @param gmt An abridged GMT object containing only the pathways that were
#' found to be significant in the ActivePathways analysis.
#' @param cytoscape_file_tag The user-defined file prefix and/or directory defining the location of the files.
#' @param col_significance A data.table object with a column 'term_id' and a column
#' for each type of omics evidence indicating whether a term was also found to be signficiant or not
#' when considering only the genes and p-values in the corresponding column of the \code{scores} matrix.
#' If term was not found, NA's are shown in columns, otherwise the relevant lists of genes are shown.
#' @param color_palette Color palette from RColorBrewer::brewer.pal to color each
#' column in the scores matrix. If NULL grDevices::rainbow is used by default.
#' @param custom_colors A character vector of custom colors for each column in the scores matrix.
#' @param color_integrated_only A character vector of length 1 specifying the color of the "combined" pathway contribution. 
#' @import ggplot2
#'
#' @return None

prepareCytoscape <- function(terms, 
                             gmt, 
                             cytoscape_file_tag, 
                             col_significance, color_palette = NULL, custom_colors = NULL, color_integrated_only = "#FFFFF0") {
  if (!is.null(col_significance)) {
    tests <- colnames(col_significance)[3:length(colnames(col_significance))]
    tests <- substr(tests, 7, 100)
    tests <- append(tests, "combined")
    
    rows <- 1:nrow(col_significance)
    
    evidence_columns = do.call(rbind, lapply(col_significance$evidence,
                                             function(x) 0+(tests %in% x)))
    colnames(evidence_columns) = tests
    
    col_significance = cbind(col_significance[,"term_id"], evidence_columns)
    
  
    if(is.null(color_palette) & is.null(custom_colors)) {
      col_colors <- grDevices::rainbow(length(tests))
    } else if (!is.null(custom_colors)){
        if (!is.null(names(custom_colors))){
          custom_colors <- custom_colors[order(match(names(custom_colors),tests))]
      }
        custom_colors <- append(custom_colors, color_integrated_only, after = match("combined",tests))
        col_colors <- custom_colors
    } else {
      col_colors <- RColorBrewer::brewer.pal(length(tests),color_palette)
    }
    col_colors <- replace(col_colors, match("combined",tests),color_integrated_only)
    if (!is.null(names(col_colors))){
      names(col_colors)[length(col_colors)] <- "combined"
    }
    
    instruct_str <- paste('piechart:',
                          ' attributelist="', 
                          paste(tests, collapse=','),
                          '" colorlist="', 
                          paste(col_colors, collapse=','), 
                          '" showlabels=FALSE', sep='')
    col_significance[, "instruct" := instruct_str]
    
    # Writing the Files
    utils::write.table(terms, 
                file=paste0(cytoscape_file_tag, "pathways.txt"), 
                row.names=FALSE, 
                sep="\t", 
                quote=FALSE)
    utils::write.table(col_significance, 
                file=paste0(cytoscape_file_tag, "subgroups.txt"), 
                row.names=FALSE, 
                sep="\t", 
                quote=FALSE)
    write.GMT(gmt, 
              paste0(cytoscape_file_tag, "pathways.gmt"))
    
    # Making a Legend
      dummy_plot = ggplot(data.frame("tests" = factor(tests, levels = tests),
                                     "value" = 1), aes(tests, fill = tests)) +
        geom_bar() +
        scale_fill_manual(name = "Contribution", values=col_colors)

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
             filename = paste0(cytoscape_file_tag, "legend.pdf"), 
             height = legend_height, 
             width = legend_width, 
             scale = 1)
    
  } else {
    utils::write.table(terms, 
                file=paste0(cytoscape_file_tag, "pathways.txt"),
                row.names=FALSE, 
                sep="\t", 
                quote=FALSE)
    write.GMT(gmt, 
              paste0(cytoscape_file_tag, "pathways.gmt"))
  }
}
