context('columnSignificance Function')


test_that('columnSignificance agrees with testing individual columns', {
    col <- colnames(dat)[1]

    pre_res1 = ActivePathways(scores = dat, gmt = gmt, background = background, cutoff = 0.1, significant = 0.05, correction.method = 'holm', return.all = T)
    res1 <- ActivePathways:::columnSignificance(scores = dat, gmt = gmt, background = background, cutoff = 0.1, significant = 0.05, correction.method = 'holm', pvals = pre_res1$adjusted.p.val)
    res2 <- ActivePathways(dat[, col, drop=FALSE], gmt, background = background, cutoff = 0.1, significant = 0.05, correction.method='holm', geneset.filter=NULL)

    # Pathways that are significant according to columnSignificance
    comp1 <- res1$term.id[which(!is.na(res1[[paste0("Genes_",col)]]))]

    expect_true(setequal(comp1, res2$term.id))
})

test_that('Column names of result if correct', {
  pre_res = ActivePathways(scores = dat, gmt = gmt, background = background, cutoff = 0.1, significant = 0.05, correction.method = 'holm', return.all = T)
  res <- ActivePathways:::columnSignificance(scores = dat, gmt = gmt, background = background, cutoff = 0.1, significant = 0.05, correction.method = 'holm', pvals = pre_res$adjusted.p.val)
  expect_equal(colnames(res), c('term.id', 'evidence', paste0("Genes_", colnames(dat))))
})
