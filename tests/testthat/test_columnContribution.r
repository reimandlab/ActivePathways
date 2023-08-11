context('columnContribution function')


test_that('Column Contribution ratio is correct', {
    col <- colnames(dat)[1]

    res <- ActivePathways(dat, gmt, significant = 1)
    res_just_column <- ActivePathways(dat[, 1, drop = FALSE], gmt, significant = 1)
    expect_equal(res_just_column$overlap, res[[paste0("Genes_", col)]])
})
