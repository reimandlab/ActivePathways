context('columnSignificance Function')


test_that('columnSignificance agrees with testing individual columns', {
    col <- colnames(dat)[1]

    res1 <- columnSignificance(dat, gmt, background, 0.1, 0.05, 'holm')
    res2 <- activeDriverPW(dat[, col, drop=FALSE], gmt, correction.method='holm', contribution=FALSE)

    # Pathways that are significant according to columnSignificance
    comp1 <- res1$term.id[which(res1[[col]] == 1)]

    expect_true(setequal(comp1, res2$term.id))
})


test_that('All results in columnSignificance are either 0 or 1', {
    res <- columnSignificance(dat, gmt, background, 0.1, 0.05, 'holm')
    res <- as.matrix(res[, -'term.id'])
    expect_true(all(apply(res, c(1,2), function(x) identical(x, 0) || identical(x, 1))))
})


test_that('Column names of result if correct', {
    res <- columnSignificance(dat, gmt, background, 0.1, 0.05, 'holm')
    expect_equal(colnames(res), c('term.id', colnames(dat)))
})
