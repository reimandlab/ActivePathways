# Hypergeometric / Fisher's exact test
# INPUT:
# genelist:       vector of gene names being tested for enrichment
# background:     vector of all gene names
# annotations;    vector of gene names annotated to the term being tested
# RETURN:
# p.val:  a number in [0,1]
hypergeometric <- function(genelist, background, annotations) {
    genelist1 <- length(which(genelist %in% annotations))
    genelist0 <- length(genelist) - genelist1
    background1 <- length(which(background %in% annotations))
    background0 <- length(background) - background1
    counts <- matrix(data=c(genelist1, genelist0, background1, background0), nrow=2, ncol=2)
    fisher.test(counts, alternative='greater')
}

# Ordered Hypergeometric / Fisher's exact test
# The hypergeometric test is run with increasingly large numbers of genes
#   starting from the top of genelist and the lowerst pvalue is returned
# INPUT:
# genelist:     vector of gene names being tested for enrichment
# background:   vector of all gene names
# annotations:  vector of gene names annotated to the term being tested
# RETURN:
# a list with the items:
#   p.val:    a number in [0,1]
#   index:    index of gene list such that genelist[1:index] gives the lowest p value
orderedHypergeometric <- function(genelist, background, annotations) {
    f <- function(ind) {
        if (genelist[ind] %in% pathway_genes) return(1)
        hypergeometric(genelist, background, annotations)$p.val
    }
    scores <- sapply(1:length(genelist), f)
    min.score <- min(scores)
    lsit(p.val=min.score, ind=max(which(scores==min_score)))
}
