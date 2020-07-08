#' Hypergemoetric test
#'
#' Perform a hypergeometric test, also known as the Fisher's exact test, on a 2x2 contingency
#' table with the alternative hypothesis 'greater' . In this application, the test finds the
#' probability that counts[1, 1] or more genes would be found to be annotated to a term (pathway),
#' assuming the null hypothesis of genes being distributed randomly to terms. 
#'
#' @param counts A 2x2 numerical matrix representing a contingency table.
#'
#' @return a p-value of enrichment of genes in a term or pathway. 
hypergeometric <- function(counts) {
    if (any(counts < 0)) stop('counts contains negative values. Something went very wrong.')
    m <- counts[1, 1] + counts[2, 1]
    n <- counts[1, 2] + counts[2, 2]
    k <- counts[1, 1] + counts[1, 2]
    x <- counts[1, 1]
    stats::phyper(x-1, m, n, k, lower.tail=FALSE)
}


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

    # Construct the counts matrix for the first which.in[1] genes
    gl <- genelist[1:which.in[1]]
    cl <- setdiff(background, gl)
    genelist0 <- length(gl) - 1
    complement1 <- length(which(cl %in% annotations))
    complement0 <- length(cl) - complement1
    counts <- matrix(data=c(1, genelist0, complement1, complement0), nrow=2)
    scores <- hypergeometric(counts)

    if (length(which.in) == 1) return(list(p.val=scores, ind=which.in[1]))

    # Update counts and recalculate score for the rest of the indeces in which.in
    # The genes in genelist[which.in[i]:which.in[i-1]] are added to the genes
    # being tested and removed from the complement. Of these, 1 will always be
    # in annotations and the rest will not. Therefore we can just modify the
    # contingency table rather than recounting which genes are in annotations
    for (i in 2:length(which.in)) {
        diff <- which.in[i] - which.in[i-1]
        counts[1, 1] <- i
        counts[2, 1] <- counts[2, 1] + diff - 1
        counts[1, 2] <- counts[1, 2] - 1
        counts[2, 2] <- counts[2, 2] - diff + 1
        scores[i] <- hypergeometric(counts)
    }

    # Return the lowest p-value and the associated index
    min.score <- min(scores)
    
    ind = which.in[max(which(scores==min.score))]
    p.val = min.score
    
	list(p.val=p.val, ind=ind)
}
