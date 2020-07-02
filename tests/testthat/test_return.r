context("Format of the returned object and output files")


test_that("Column names of data.table is correct", {
    expect_equal(colnames(run_ap_short(dat)),
                 c('term.id', 'term.name', 'adjusted.p.val', 'term.size', 'overlap'))
    expect_equal(colnames(run_ap_short_contribution(dat)),
                 c('term.id', 'term.name', 'adjusted.p.val', 'term.size', 'overlap', 
                 'evidence', 'Genes_cds', 'Genes_promoter', 'Genes_enhancer'))
})


test_that("All results or only significant ones are returned", {
    expect_true(all(ActivePathways(dat, gmt, cytoscape.file.tag = NA)$p.val < 0.05))
    expect_equal(nrow(ActivePathways(dat, gmt, cytoscape.file.tag = NA, significant = 1)), length(gmt))
})


test_that("No significant results are found", {
    expect_warning(res1 <- ActivePathways(dat, gmt, cytoscape.file.tag = NA, significant = 0),
                   "No significant terms were found.", fixed = TRUE)
    expect_equal(nrow(res1), NULL)
})


