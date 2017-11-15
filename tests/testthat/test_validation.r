context("Validation on parameters")

# Prepare testing data
library(mpea)
gmt <- read.GMT('test.gmt')
dat <- as.matrix(read.table('test_data.txt', header=TRUE, row.names='Gene'))
dat[is.na(dat)] <- 1

file.names <- paste(c('terms', 'groups', 'smallgmt'), '.txt', sep="")

test_that("cytoscape.filenames specified", {
    expect_error(mpea(dat, gmt, contribution=TRUE, cytoscape.filenames=file.names[1]),
                 "Must supply 3 file names to cytoscape.filenames", fixed=TRUE)
    expect_error(mpea(dat, gmt, contribution=FALSE, cytoscape.filenames=file.names[1]),
                 "Must supply 2 file names to cytoscape.filenames", fixed=TRUE)
    expect_error(mpea(dat, gmt, contribution=FALSE, significant=1, cytoscape.filenames=file.names[c(1,3)]), NA)
    expect_message(mpea(dat, gmt, contribution=FALSE, significant=1, cytoscape.filenames=file.names),
                   "Column contributions will not be evaluated so the contribution matrix is not being written. cytoscape.filenames[2] will be ignored", fixed=TRUE)
})

test_that("significant is valid", {
    expect_error(mpea(dat, gmt, significant=-0.1), "significant must be a value in [0,1]", fixed=TRUE)
    expect_error(mpea(dat, gmt, significant = 1.1), "significant must be a value in [0,1]", fixed=TRUE)
    expect_warning(mpea(dat, gmt, significant=0, return.all=TRUE), "No significant terms were found")
    expect_error(mpea(dat, gmt, significant=1, return.all=TRUE), NA)
})

test_that("cutoff is valid", {
    expect_error(mpea(dat, gmt, cutoff=-0.1), "cutoff must be a value in [0,1]", fixed=TRUE)
    expect_error(mpea(dat, gmt, cutoff = 1.1), "cutoff must be a value in [0,1]", fixed=TRUE)
    expect_error(mpea(dat, gmt, cutoff=0, return.all=TRUE), "No genes made the cutoff", fixed=TRUE)
    expect_error(mpea(dat, gmt, cutoff=1, return.all=TRUE), NA)
})

test_that("genes not found in background are removed", {
    expect_message(mpea(dat, gmt, background=rownames(dat)[1], significant=1, cutoff=1), "99 rows were removed from scores because they are not found in the background")
    expect_error(mpea(dat, gmt, background='qwerty'), "scores does not contain any genes in the background")
})
