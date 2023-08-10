#' Ordered Hypergeometric Test
#'
#' Perform a series of hypergeometric tests (a.k.a. Fisher's Exact tests), on a ranked list of genes ordered
#' by significance against a list of annotation genes. The hypergeometric tests are executed with 
#' increasingly larger numbers of genes representing the top genes in order of decreasing scores. 
#' The lowest p-value of the series is returned as the optimal enriched intersection of the ranked list of genes
#' and the biological term (pathway). 
#'
#' @param genelist Character vector of gene names, assumed to be ordered by decreasing importance. 
#' For example, the genes could be ranked by decreasing significance of differential expression. 
#' @param background Character vector of gene names. List of all genes used as a statistical background (i.e., the universe)
#' @param annotations Character vector of gene names. A gene set representing a functional term, process or biological pathway. 
#'
#' @return a list with the items:
#'   \describe{
#'     \item{p.val}{The lowest obtained p-value}
#'     \item{ind}{The index of \code{genelist} such that \code{genelist[1:ind]}
#'       gives the lowest p-value}
#'  }
#' @export
#'
#' @examples
#'    orderedHypergeometric(c('HERC2', 'SP100'), c('PHC2', 'BLM', 'XPC', 'SMC3', 'HERC2', 'SP100'),
#'                          c('HERC2', 'PHC2', 'BLM'))
orderedHypergeometric <- function(genelist, background, annotations) {
  # Only test subsets of genelist that end with a gene in annotations since
  # these are the only tests for which the p-value can decrease
  which.in <- which(genelist %in% annotations)
  if (length(which.in) == 0) return(list(p.val=1, ind=1))
  
  scores <- c()
  
  #go through all the lists and calculate the hyper geometric test
  for (i in 1:length(which.in)) {
    scores[i] <- 
      stats::phyper(
        q = i - 1, #overlap - 1 
        m = length(annotations),
        n = length(background) - length(annotations),
        k = which.in[i],
        lower.tail = FALSE)
  }
  
  # Return the lowest p-value and the associated index
  min.score <- min(scores)
  
  ind = which.in[max(which(scores==min.score))]
  p.val = min.score
  
  list(p.val=p.val, ind=ind)
}