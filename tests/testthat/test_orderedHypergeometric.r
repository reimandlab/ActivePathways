context("Ordered Hypergeometric Statistical Test")


test_that('hypergeometric gives the same results as fisher.test', {
    counts <- matrix(c(0,0,0,0), nrow=2)
    expect_equal(hypergeometric(counts), fisher.test(counts, alternative='greater')$p.value)

    counts <- matrix(c(0, 5, 16, 1683), nrow=2)
    expect_equal(hypergeometric(counts), fisher.test(counts, alternative='greater')$p.value)

    counts <- matrix(c(2, 2, 2, 2), nrow=2)
    expect_equal(hypergeometric(counts), fisher.test(counts, alternative='greater')$p.value)
})


test_that('orderedHypergeometric returns the lowest p-value and correct index', {
    genelist <- c('HERC2', 'SP100', 'BLM')
    background <- c('PHC2', 'BLM', 'XPC', 'SMC3', 'HERC2', 'SP100')
    annotations <- c('HERC2', 'PHC2', 'BLM')

    get_pvalue <- function(genes) {
        complement <- setdiff(background, genes)
        genelist1 <- length(which(genes %in% annotations))
        genelist0 <- length(genes) - genelist1
        complement1 <- length(which(complement %in% annotations))
        complement0 <- length(complement) - complement1
        counts <- matrix(c(genelist1, genelist0, complement1, complement0), 2)
        hypergeometric(counts)
    }
    p.values <- sapply(1:length(genelist), function(i) get_pvalue(genelist[1:i]))

    smallest.value <- min(p.values)
    smallest.index <- match(smallest.value, p.values)
    exp <- list(p.val=smallest.value, ind=smallest.index)

    expect_equal(orderedHypergeometric(genelist, background, annotations), exp)
})
