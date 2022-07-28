context("Test the enrichmentAnalysis function")


test_that('Overlap Found by enrichmentAnalysis is correct', {
    res <- enrichmentAnalysis(ea_genelist, ea_gmt, ea_background)
    res2 <- enrichmentAnalysis(ea_genelist2, ea_gmt, ea_background2)

    expect_true(setequal(res[[1, 'overlap']], c('BLM')))
    expect_true(setequal(res[[2, 'overlap']], c('HERC2', 'SP100', 'BLM')))
    expect_true(setequal(res2[[3, 'overlap']], c('HERC2')))
    expect_true(setequal(res[[4, 'overlap']], NULL))
})
