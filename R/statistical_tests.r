#' Hypergeometric Test
#'
#' Perform a Hypergeometric, aka Fisher's exact test, on a list of significant
#' genes and annotation genes
#'
#' @param complement character vector of gene names. All other genes not being
#'   tested for enrichment. Ie \code{setdiff(background, genelist)}
#'
#' @inheritParams orderedHypergeometric
#'
#' @return a p-value
#'
#' @examples
#' \dontrun{
#'   hypergeometric(c('HERC2', 'SP100'), c('PHC2', 'BLM', 'NDC1', 'XPC', 'SMC3'), c('HERC2', 'PHC2', 'BLM'))
#' }
hypergeometric <- function(genelist, complement, annotations) {
    genelist1 <- length(which(genelist %in% annotations))
    genelist0 <- length(genelist) - genelist1
    complement1 <- length(which(complement %in% annotations))
    complement0 <- length(complement) - complement1
    counts <- matrix(data=c(genelist1, genelist0, complement1, complement0), nrow=2)
    fisher.test(counts, alternative='greater')
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
#'  @examples
#'  \dontrun{
#'    orderedHypergeometric(c('HERC2', 'SP100'), c('PHC2', 'BLM', 'NDC1', 'XPC', 'SMC3', 'HERC2', 'SP100'),
#'                          c('HERC2', 'PHC2', 'BLM'))
#'  }
orderedHypergeometric <- function(genelist, background, annotations) {
    # remove genes in genelist from background to form a list of all other genes
    # genes that are not being tested in each call will be added back
    complement <- setdiff(background, genelist)

    # Test only when a gene is in annotation as that's the only case that
    # could reduce the p-value
    which.in <- which(genelist %in% annotations)

    if (length(which.in) == 0) return(list(p.val=1, ind=1))

    # Test first i genes in genelist. The rest of genelist has to be added to
    # complement
    f <- function(i) {
        hypergeometric(genelist[1:i],
                       c(genelist[(i+1):length(genelist)], complement),
                       annotations)
    }
    scores <- sapply(which.in, function(i) f(i)$p.val)

    min.score <- min(scores)
    list(p.val=min.score, ind=which.in[max(which(scores==min.score))])
}
