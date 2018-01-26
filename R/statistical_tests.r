#' Hypergemoetric Test
#'
#' Perform a Hypergeometric test, aka Fisher's exact test, on a 2x2 contingency
#' table with a 'greater' alternative hypothesis. That is, find the
#' probability that counts[1, 1] or more genes would be found in annotations,
#' assuming the null hypothesis.
#'
#' @param counts a 2x2 numerical matrix representing a contingency table
#'
#' @return a p-value
#'
#' @examples
#' \dontrun{
#'   hypergeometric(matrix(data=c(1, 8, 14, 836), nrow=2))
#' }
hypergeometric <- function(counts) {
    if (any(counts < 0)) stop('counts contains negative values. Something went very wrong.')
    m <- counts[1, 1] + counts[2, 1]
    n <- counts[1, 2] + counts[2, 2]
    k <- counts[1, 1] + counts[1, 2]
    x <- counts[1, 1]
    phyper(x-1, m, n, k, lower.tail=FALSE)
}


#' Ordered Hypergeometric Test
#'
#' Perform a Hypergeometric, aka Fisher's Exact test, on a list of genes ordered
#' by significance against a list of annotation genes
#'
#' The hypergeometric test is run with increasingly large numbers of genes
#' starting from the top, and the lowest p-value is returned
#'
#' @param genelist character vector of gene names. List of differentially
#'   expressed genes being tested for enrichment
#' @param background character vector of gene names. List of all genes being
#'   used as a statistical background
#' @param annotations character vector of gene names. List of genes annotated to
#'   the term being tested.
#'
#' @return a list with the items:
#'   \describe{
#'     \item{p.val}{The lowest obtained p-value}
#'     \item{ind}{The index of \code{genelist} such that \code{genelist[1:ind]}
#'       gives the lowest p-value}
#'  }
#'
#' @examples
#' \dontrun{
#'    orderedHypergeometric(c('HERC2', 'SP100'), c('PHC2', 'BLM', 'XPC', 'SMC3', 'HERC2', 'SP100'),
#'                          c('HERC2', 'PHC2', 'BLM'))
#' }
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

    list(p.val=min.score*length(scores), ind=which.in[max(which(scores==min.score))])
#    list(p.val=min.score, ind=which.in[max(which(scores==min.score))])
}
