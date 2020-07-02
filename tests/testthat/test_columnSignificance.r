context('columnSignificance Function')


test_that('columnSignificance agrees with testing individual columns', {
    col <- colnames(dat)[1]

    res1 <- columnSignificance(dat, gmt, background, 0.1, 0.05, 'holm', rep(0.05, length(gmt)))
    res2 <- ActivePathways(dat[, col, drop = FALSE], gmt, 
    		correction.method = 'holm', geneset.filter = NULL)

    # Pathways that are significant according to columnSignificance
	comp1 = res1$term.id[ sapply(res1[["Genes_cds"]], function(x) !all(is.na(x))) ]

    expect_true(setequal(comp1, res2$term.id))
})

test_that('Column names of columnSignificance result is correct', {
    res <- columnSignificance(dat, gmt, background, 0.1, 0.05, 'holm', rep(0.05, length(gmt)))
    expect_equal(colnames(res), c('term.id', 'evidence', paste0("Genes_", colnames(dat))))
})
