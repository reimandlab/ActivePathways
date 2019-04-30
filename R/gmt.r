#' Read and Write GMT files
#'
#' Functions to read and Write Gene Matrix Transposed (GMT) files and to test if
#' an object inherits from GMT
#'
#' A GMT file describes gene sets, such as pathways. GMT files are tab delimeted
#' and each row contains a term id, a term name, and all genes annotated to the
#' term.
#'
#' @format
#' A GMT object is a named list of terms, where each term is a list with the items:
#' \describe{
#'     \item{id}{The term id}
#'     \item{name}{The full name of the term}
#'     \item{genes}{A character vector of genes annotated to this term}
#'   }
#' @exportClass GMT
#' @rdname GMT
#' @name GMT
#' @aliases GMT gmt
#'
#' @param filename Location of the gmt file
#' @param gmt a GMT object
#' @param x object to test
#' @param i index of GMT object
#'
#' @return \code{read.GMT} returns a GMT object. \cr
#' \code{write.GMT} returns NULL. \cr
#' \code{is.GMT} returns TRUE if \code{x} is a GMT object, else FALSE
#'
#'
#' @examples
#' gmt <- read.GMT(system.file('extdata', 'hsapiens_REAC_subset.gmt', package='ActivePathways'))
#' is.GMT(gmt)
#' gmt[1:10]
#' gmt[[1]]
#' gmt[1]$id
#' gmt[1]$genes
#' gmt[1]$name
#' gmt$`REAC:3108214`
#' \donttest{
#'   write.GMT(gmt, 'filename.gmt')
#' }
NULL

#' Reading GMT files
#' @rdname GMT
#' @export
read.GMT <- function(filename) {
    gmt <- strsplit(readLines(filename), '\t')
    names(gmt) <- sapply(gmt, `[`, 1)
    gmt <- lapply(gmt, function(x) { list(id=x[1], name=x[2], genes=x[-c(1,2)]) })
    class(gmt) <- 'GMT'
    gmt
}

#' Writing GMT files
#' @rdname GMT
#' @export
write.GMT <- function(gmt, filename) {
    if (!is.GMT(gmt)) stop("gmt is not a valid GMT object")
    sink(filename)
    for (term in gmt) {
        cat(term$id, term$name, paste(term$genes, collapse="\t"), "\n", sep="\t")
    }
    sink()
}

#####  Subsetting functions #####
# Treat as a list but return an object of "GMT" class
#' @rdname GMT
#' @export
`[.GMT` <- function(x, i) {
    x <- unclass(x)
    res <- x[i]
    class(res) <- c('GMT')
    res
}
#' @rdname GMT
#' @export
`[[.GMT` <- function(x, i) {
    x <- unclass(x)
    x[[i, exact = TRUE]]
}
#' @rdname GMT
#' @export
`$.GMT` <- function(x, i) {
    x[[i]]
}

#' Checks whether an object is a GMT
#' @rdname GMT
#' @export
is.GMT <- function(x) inherits(x, 'GMT')

#' Make a background list of genes
#'
#' Returns a character vector of all genes in a GMT object
#'
#' @param gmt a \link{GMT} object
#' @return a character vector containing all genes in GMT
#' @export
#'
#' @examples
#' gmt <- read.GMT(system.file('extdata', 'hsapiens_REAC_subset.gmt', package='ActivePathways'))
#' makeBackground(gmt)
makeBackground <- function(gmt) {
  if (!is.GMT(gmt)) stop('gmt is not a valid GMT object')
  unlist(Reduce(function(x, y) union(x, y$genes), gmt, gmt[[1]]$genes))
}
