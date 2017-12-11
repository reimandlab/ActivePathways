context("Validation of Cytoscape Files and Test that the files are written")


test_that("cytoscape.filenames specified", {
    expect_error(activeDriverPW(dat, gmt, contribution=TRUE, cytoscape.filenames=file.names[1]),
                 "Must supply 3 file names to cytoscape.filenames", fixed=TRUE)
    expect_error(activeDriverPW(dat, gmt, contribution=FALSE, cytoscape.filenames=file.names[1]),
                 "Must supply 2 file names to cytoscape.filenames", fixed=TRUE)
    expect_error(activeDriverPW(dat, gmt, contribution=FALSE, significant=1,
                                cytoscape.filenames=file.names[c(1,3)]), NA)
    expect_message(activeDriverPW(dat, gmt, contribution=FALSE, significant=1,
                                  cytoscape.filenames=file.names),
                   "Column contributions will not be evaluated so the contribution matrix is not being written. cytoscape.filenames[2] will be ignored", fixed=TRUE)
})


test_that("Cytoscape files are written", {
    suppressWarnings(file.remove(file.names))
    activeDriverPW(dat, gmt, cytoscape.filenames=file.names, contribution=TRUE,
                   significant=0.9, cutoff=1)
    expect_equal(file.exists(file.names), c(TRUE, TRUE, TRUE))

    suppressWarnings(file.remove(file.names))
    activeDriverPW(dat, gmt, cytoscape.filenames=file.names, contribution=FALSE,
                   significant=0.9, cutoff=1)
    expect_equal(file.exists(file.names), c(TRUE, FALSE, TRUE))

    suppressWarnings(file.remove(file.names))
    suppressWarnings(activeDriverPW(dat, gmt, cytoscape.filenames=file.names,
                                    contribution=TRUE, significant=0, return.all=TRUE))
    expect_equal(file.exists(file.names), c(FALSE, FALSE, FALSE))

    suppressWarnings(file.remove(file.names))
    suppressWarnings(activeDriverPW(dat, gmt, cytoscape.filenames=NULL, contribution=TRUE))
    expect_equal(file.exists(file.names), c(FALSE, FALSE, FALSE))

    suppressWarnings(file.remove(file.names))
})
