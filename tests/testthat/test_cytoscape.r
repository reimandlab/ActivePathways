context("Validation of Cytoscape Files and Test that the files are written")


test_that("cytoscape.filenames specified", {
    expect_error(ActivePathways(dat, gmt, cytoscape.file.dir=1),
                 "cytoscape.file.dir must be a string", fixed=TRUE)
    expect_error(ActivePathways(dat, gmt, significant=1,
                                cytoscape.file.dir=tempdir()), NA)
})


test_that("Cytoscape files are written", {
    suppressWarnings(file.remove(file.names))
    ActivePathways(dat, gmt, cytoscape.file.dir=tempdir(),
                   significant=0.9, cutoff=1)
    expect_equal(file.exists(file.names), c(TRUE, TRUE, TRUE, TRUE))

    suppressWarnings(file.remove(file.names))
    suppressWarnings(ActivePathways(dat, gmt, cytoscape.file.dir=tempdir(),
                                    significant=0, return.all=TRUE))
    expect_equal(file.exists(file.names), c(FALSE, FALSE, FALSE, FALSE))

    suppressWarnings(file.remove(file.names))
    suppressWarnings(ActivePathways(dat, gmt, cytoscape.file.dir=NULL))
    expect_equal(file.exists(file.names), c(FALSE, FALSE, FALSE, FALSE))

    suppressWarnings(file.remove(file.names))
})
