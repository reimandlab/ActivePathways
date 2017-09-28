# Read in a gmt file
# 
# INPUT: 
# filename: Location of the gmt file
#
# RETURN:
# A gmt object, which is a list of term ids where each termid is a list with the items:
#   id:     term id
#   name:   full name of the term
#   genes:  character vector of genes annotated to this term

read_gmt <- function(filename) {
    gmt <- strsplit(readLines(filename), '\t')
    names(gmt) <- sapply(gmt, `[`, 1)
    gmt <- lapply(gmt, function(x) { list(id=x[1], name=x[2], genes=x[-c(1,2)]) })
    class(gmt) <- 'gmt'
    gmt
}


# Create a background list of all genes in the gmt object
# INPUT:
# gmt:      a gmt object
# RETURN: 
# a character vector of gene names
make_background <- function(gmt) {
    unlist(Reduce(function(x, y) union(x, y$genes), gmt, gmt[[1]]$genes))
}


# Pretty-print the gmt file
print.gmt(gmt) <- function(x) {
    for (term in gmt) {
        cat(paste(term$id, "-", term$name, "\n"))
        if (length(term$genes) > 10) {
            cat(term$genes[1:10])
            cat("...")
        } else {
            cat(term$genes)
        }
        cat("\n\n")
    }
}
