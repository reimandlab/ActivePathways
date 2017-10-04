#' Read a gmt file
#'
#' Returns a gmt object read in from a file
#'
#' @param filename a character vector containing the location of the gmt file
#' @return a gmt object, which is a lest of terms where each term is a list with the items:
#'   \itemize{
#'     \item{"id"}{The term id}
#'     \item{"name"}{The full name of the term}
#'     \item{"genes"}{A character vector of genes annotated to this term}
#'   }
#' @exportClass gmt
#' @export
readGMT <- function(filename) {
    gmt <- strsplit(readLines(filename), '\t')
    names(gmt) <- sapply(gmt, `[`, 1)
    gmt <- lapply(gmt, function(x) { list(id=x[1], name=x[2], genes=x[-c(1,2)]) })
    class(gmt) <- 'GMT'
    gmt
}

#' Make a background list of genes
#'
#' Returns a character vector of all genes in a gmt object
#'
#' @param gmt a gmt object (see \code{\link{readGMT}})
#' @return a character vector containing all genes in gmt
#' @export
makeBackground <- function(gmt) {
    if (!is.GMT(gmt)) stop('gmt is not a valid GMT object')
    unlist(Reduce(function(x, y) union(x, y$genes), gmt, gmt[[1]]$genes))
}

### Subsetting functions. Treat as a list but return an object of "GMT" class
`[.GMT` <- function(x, i) {
    x <- unclass(x)
    res <- x[i]
    class(res) <- c('GMT')
    res
}
`[[.GMT` <- function(x, i, exact = TRUE) {
    x <- unclass(x)
    x[[i, exact = exact]]
}
`$.GMT` <- function(x, i) {
    x[[i]]
}
is.GMT <- function(x) inherits(x, 'GMT')

#' Print a gmt object
#'
#' @export
print.GMT <- function(x, ...) {
    num.lines <- min(length(x), getOption("max.print", 99999))
    num.trunc <- length(x) - num.lines
    str <- sapply(x[1:num.lines], function(a) paste(a$id, "-", a$name, "\n", paste(a$genes, collapse=", "), '\n\n'))
    cat(str)
    if (num.trunc == 1) {
        cat('[ reached getOption("max.print") -- omitted 1 term ]')
    } else if (num.trunc > 1) {
        cat(paste('[ reached getOption("max.print") -- ommitted', num.trunc, 'terms ]'))
    }
}
