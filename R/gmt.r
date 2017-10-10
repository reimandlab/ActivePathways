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
#'
#' @return \code{read.GMT} returns a GMT object. \cr
#' \code{write.GMT} returns NULL. \cr
#' \code{is.GMT} returns TRUE if \code{x} is a GMT object, else FALSE
#'
NULL

#' @rdname GMT
#' @export
read.GMT <- function(filename) {
    gmt <- strsplit(readLines(filename), '\t')
    names(gmt) <- sapply(gmt, `[`, 1)
    gmt <- lapply(gmt, function(x) { list(id=x[1], name=x[2], genes=x[-c(1,2)]) })
    class(gmt) <- 'GMT'
    gmt
}

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

#' Make a background list of genes
#'
#' Returns a character vector of all genes in a GMT object
#'
#' @param gmt a \link{GMT} object
#' @return a character vector containing all genes in GMT
#' @export
makeBackground <- function(gmt) {
    if (!is.GMT(gmt)) stop('gmt is not a valid GMT object')
    unlist(Reduce(function(x, y) union(x, y$genes), gmt, gmt[[1]]$genes))
}

### Subsetting functions. Treat as a list but return an object of "GMT" class
#' @export
`[.GMT` <- function(x, i) {
    x <- unclass(x)
    res <- x[i]
    class(res) <- c('GMT')
    res
}
#' @export
`[[.GMT` <- function(x, i, exact = TRUE) {
    x <- unclass(x)
    x[[i, exact = exact]]
}

#' @export
`$.GMT` <- function(x, i) {
    x[[i]]
}

#' @export
#' @rdname GMT
is.GMT <- function(x) inherits(x, 'GMT')

# Print a GMT object
#' @export
print.GMT <- function(x, ...) {
    num.lines <- min(length(x), getOption("max.print", 99999))
    num.trunc <- length(x) - num.lines
    cat(sapply(x[1:num.lines], function(a) paste(a$id, "-", a$name, "\n",
                                                 paste(a$genes, collapse=", "), '\n\n')))
    if (num.trunc == 1) {
        cat('[ reached getOption("max.print") -- omitted 1 term ]')
    } else if (num.trunc > 1) {
        cat(paste('[ reached getOption("max.print") -- ommitted', num.trunc, 'terms ]'))
    }
}
