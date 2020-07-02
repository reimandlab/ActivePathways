context("Validation of Cytoscape Files and Test that the files are written")


test_that("cytoscape.filenames specified", {

    expect_error(ActivePathways(dat, gmt, cytoscape.file.tag = CStag, significant = 1), NA)
    expect_message(ActivePathways(dat[,1, drop = F], gmt, cytoscape.file.tag = CStag, significant = 1),
			"Scores matrix contains only one column. Column contributions will not be calculated.", 
			fixed = TRUE)
})


test_that("Cytoscape files are written", {
	
	CS_fnames = paste0(CStag, c("pathways.txt", "subgroups.txt", "pathways.gmt", "legend.pdf"))
	
    suppressWarnings(file.remove(CS_fnames))
    ActivePathways(dat, gmt, cytoscape.file.tag = CStag, significant = 0.9, cutoff = 1)
    expect_equal(file.exists(CS_fnames), c(TRUE, TRUE, TRUE, TRUE))

    suppressWarnings(file.remove(CS_fnames))
    ActivePathways(dat[,1, drop = F], gmt, cytoscape.file.tag = CStag, significant = 0.9, cutoff = 1)
    expect_equal(file.exists(CS_fnames), c(TRUE, FALSE, TRUE, FALSE))

    suppressWarnings(file.remove(CS_fnames))
 	suppressWarnings(ActivePathways(dat, gmt, cytoscape.file.tag = CStag, significant = 0))
    expect_equal(file.exists(CS_fnames), c(FALSE, FALSE, FALSE, FALSE))

    suppressWarnings(file.remove(CS_fnames))
    suppressWarnings(ActivePathways(dat, gmt, cytoscape.file.tag = NA))
    expect_equal(file.exists(CS_fnames), c(FALSE, FALSE, FALSE, FALSE))

    suppressWarnings(file.remove(CS_fnames))
})
