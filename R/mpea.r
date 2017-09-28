# INPUT:
# scores:       a numerical matrix where each row is a gene and each column is a test.
#               rownames must be the genes. All values should be 0<=p<=1 and NAs should be converted to 1
# gmt:          gmt object to be used for enrichment analysis. If a filename, the object will be created from file
# background:   vector of gene names to use as a statistical background. By default, all genes in the gmt object
# cutoff:       maximum p-value of genes to be used for enrichment analysis. Each gene in scores will receive a single 
#               p-value by merging the p-values from all tests. The genes will then be ordered by p-value and any genes
#               with a p-value greater than cutoff will be discarded
# significant:  If a number in [0, 1], only terms with a p-value smaller than significant will be returned.
#               If NA, FALSE, etc, all terms will be returned
# return.all:   Whether to return results for all terms. If FALSE, only terms with p.val < significant will be returned
# cytoscape.filenames:    a vector with 3 filenames to write data for cytoscape to. If FALSE, do not write these files.
#       [1]     List of significant terms
#       [2]     Matrix of column contributions to each term
#       [3]     Shortened version of the gmt file, containing only terms 
mpea <- function(scores, gmt, merge.method=c("Fisher", "Browns"), correction.method=c("fdr"), 
                 background=makeBackground(gmt), cutoff=0.1, significant=0.05, return.all=FALSE,
                 cytoscape.filenames=c('termlist.txt', 'subgroups.txt', 'abridged.gmt')) {
    # Type-checking
    merge.method <- match.arg(merge.method)
    correction.method <- match.arg(correction.method)
    if (!is.numeric(scores)) stop("scores must be a numeric matrix")
    if (!is.numeric(cutoff) || cutoff < 0 || cutoff > 1) stop("cutoff must be in [0,1]")
    if (!is.numeric(significant) || significant < 0 || significant > 0) { 
        stop("significant must be a logical value or in [0,1]")
    }
    if (!is.character(background)) stop("background must be a character vector")
    if (class(gmt) != 'gmt') gmt <- read_gmt(gmt)
    

    merged.scores <- merge_p_values(scores)
    merged.scores <- merged.scores[merged.scores <= cutoff]
    ordered.scores <- names(merged.scores)[order(merged.scores)]
    
    res <- gsea(ordered.scores, gmt, background)
    contribution <- columnContribution(scores, gmt, background, significant, cutoff, correction.method)
    res <- cbind(res, contribution[, -'term.name'])

    significant.indeces <- which(res$p.val <= significant)

    if (cytoscape.filenames) prepareCytoscape(contribution[significant.indeces], gmt, cytoscape.filenames)

    if (return.all) return(res)
    return(res[significant.indexes])
}


# Perform gsea on an ordered list of genes
gsea <- function(genelist, gmt, background) {
    dt <- data.table(term.id=names(gmt))
    
    for (i in 1:length(gmt)) {
        term <- gmt[[i]]
        tmp <- orderedHypergeometric(term$genes, scores, background)
        overlap <- scores[scores[1:tmp$ind] %in% term$genes]
        if (length(overlap) == 0) overlap <- NA
        set(dt, i, 'term.name', term$name)
        set(dt, i, 'p.val', tmp$p.val)
        set(dt, i, 'term.size', length(term$genes))
        set(dt, i, 'query.size', length(scores))
        set(dt, i, 'overlap', list(list(overlap)))
    }
    dt
}

# Get contribution of each column to each term
columnContribution <- function(scores, gmt, background, significant, cutoff, correction.method) {
    dt <- data.table(term.id=names(gmt))

    for (col in colnames(scores)) {
        col.scores <- p.adjust(scores[, col, drop=TRUE], method=correction.method)
        col.scores <- col.scores[col.scores <= cutoff]
        ordered.scores <- scores[order(scores[, col]), col]

        p.vals <- sapply(term.id, function(x) orderedHypergeometric(gmt[[x]]$genes, ordered.scores, background)$p.val)
        p.vals <- sapply(p.vals, function(x) ifelse(x<=significant, 1, 0))
        dt[, (col) := p.vals]
     }
    dt
}

# Prepare files for building an Enrichment Map in Cytoscape
# INPUT:
# enrichment.matrix:    a logical matrix where each row is a term and each column is a test
#                       a TRUE indicates that a pathway is enriched considering that test only
# gmt:                  gmt object
# filenames:            a vector with 3 filenames to write data for cytoscape to
#       [1]             List of significant terms
#       [2]             Matrix of column contributions to each term
#       [3]             Shortened version of the gmt file, containing only terms 
prepareCytoscape <- function(enrichment.matrix, filenames) {
    termlist <- enrichment_matrix[, .(term.id, term.name)]
    termlist[, p.val := 0.05]

    subgroups <- enrichment.matrix[, -'term.name']
    tests <- colnames(subgroups[, -'term.id'])
    instruct.str <- paste('pichart:', 
                          ' attributelist="', paste(tests,collapse=","), '"',
                          ' colorlist="', paste(rainbow(length(tests)), collapse=","), '"',
                          ' showlabels=FALSE', sep="")
    subgroups[, instruct := instruct.str]
    
    abridged.gmt <- gmt[enrichment.matrix[, 'term.id']]

    write.table(termlist, file=filenames[1], row.names=FALSE, sep="\t", quote=FALSE)
    write.table(subgroups, file=filenames[2], row.names=FALSE, sep="\t", quote=FALSE)
    write.table(abridged.gmt, file=filenames[3], row.names=FALSE, sep="\t", quote=FALSE)
}
    


