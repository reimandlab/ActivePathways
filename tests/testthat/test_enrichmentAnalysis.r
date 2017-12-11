context("Test the enrichmentAnalysis function")


test_that('Overlap Found by enrichmentAnalysis is correct', {
    res <- enrichmentAnalysis(ea.genelist, ea.gmt, ea.background)
    res2 <- enrichmentAnalysis(ea.genelist2, ea.gmt, ea.background2)

    expect_true(setequal(res[[1, 'overlap']], c('BLM')))
    expect_true(setequal(res[[2, 'overlap']], c('HERC2', 'SP100', 'BLM')))
    expect_true(setequal(res2[[3, 'overlap']], c('HERC2')))
    expect_true(setequal(res[[4, 'overlap']], NA))
})
