context('columnContribution function')


test_that('Column Contribution ratio is correct', {
    col <- colnames(dat)[1]

    res <- activePathways(dat, gmt, return.all=TRUE)
    res.without.column <- activePathways(dat[, -1, drop=FALSE], gmt, return.all=TRUE)
    expected.contribution <- -log10(res$p.val / res.without.column$p.val)

    expect_equal(res[[col]], expected.contribution)
})
