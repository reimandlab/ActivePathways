context("Validation on input to ActivePathways")


test_that("scores is a numeric matrix with valid p-values", {
  dat2 <- dat
  dat2[1, 1] <- 'a'
  expect_error(run_ap_short(dat2), 'scores must be a numeric matrix')

  dat2 <- dat
  dat2[1, 1] <- NA
  expect_error(run_ap_short(dat2), 'scores may not contain missing values')

  dat2[1, 1] <- -0.1
  expect_error(run_ap_short(dat2), "All values in scores must be in [0,1]", fixed=TRUE)

  dat2[1, 1] <- 1.1
  expect_error(run_ap_short(dat2), "All values in scores must be in [0,1]", fixed=TRUE)

  dat2[1, 1] <- 1
  expect_error(run_ap_short(dat2), NA)

  dat2[1, 1] <- 0
  expect_error(run_ap_short(dat2), NA)
})

test_that("significant is valid", {
    expect_error(ActivePathways(dat, gmt, significant=-0.1),
                 "significant must be a value in [0,1]", fixed=TRUE)
    expect_error(ActivePathways(dat, gmt, significant = 1.1),
                 "significant must be a value in [0,1]", fixed=TRUE)
    expect_error(ActivePathways(dat, gmt, significant=NULL),
                 "length(significant) == 1 is not TRUE", fixed=TRUE)
    expect_error(ActivePathways(dat, gmt, significant=c(1,2)),
                 "length(significant) == 1 is not TRUE", fixed=TRUE)
    expect_error(ActivePathways(dat, gmt, significant='qwe'),
                 "is.numeric(significant) is not TRUE", fixed=TRUE)
    expect_warning(ActivePathways(dat, gmt, significant = 0),
                   "No significant terms were found")
    expect_error(ActivePathways(dat, gmt, significant = 1), NA)
})


test_that("cutoff is valid", {
    expect_error(ActivePathways(dat, gmt, cutoff=-0.1),
                 "cutoff must be a value in [0,1]", fixed=TRUE)
    expect_error(ActivePathways(dat, gmt, cutoff = 1.1),
                 "cutoff must be a value in [0,1]", fixed=TRUE)
    expect_error(ActivePathways(dat, gmt, cutoff=NULL),
                 "length(cutoff) == 1 is not TRUE", fixed=TRUE)
    expect_error(ActivePathways(dat, gmt, cutoff=c(1,2)),
                 "length(cutoff) == 1 is not TRUE", fixed=TRUE)
    expect_error(ActivePathways(dat, gmt, cutoff='qwe'),
                 "is.numeric(cutoff) is not TRUE", fixed=TRUE)
    expect_error(ActivePathways(dat, gmt, cutoff=0),
                 "No genes made the cutoff", fixed=TRUE)
    expect_error(ActivePathways(dat, gmt, cutoff=1), NA)
})


test_that("background is a character vector", {
    error_msg <- "background must be a character vector"
    expect_error(ActivePathways(dat, gmt, background=c(1,5,2)), error_msg)
    expect_error(ActivePathways(dat, gmt, background=matrix(c('a', 'b', 'c', 'd'), 2)), error_msg)
})


test_that("genes not found in background are removed", {
    expect_message(ActivePathways(dat, gmt, background=rownames(dat)[-(1:10)], significant=1, cutoff=1),
                   "10 rows were removed from scores because they are not found in the background")
    expect_error(ActivePathways(dat, gmt, background='qwerty'),
                 "scores does not contain any genes in the background")
})

test_that("geneset.filter is a numeric vector of length 2", {
    expect_error(ActivePathways(dat, gmt, geneset.filter=1), 
                 "geneset.filter must be length 2")
    expect_error(ActivePathways(dat, gmt, geneset.filter=list(1,2)),
                 "geneset.filter must be a numeric vector")
    expect_error(ActivePathways(dat, gmt, geneset.filter=c('q', 2)),
                 "geneset.filter must be a numeric vector")
    expect_error(ActivePathways(dat, gmt, geneset.filter=c(1, -2)),
                 "geneset.filter limits must be positive")
    expect_error(ActivePathways(dat, gmt, geneset.filter=c(0, 0)), 
                 "No pathways in gmt made the geneset.filter", fixed=TRUE)
    expect_message(ActivePathways(dat, gmt, geneset.filter=c(NA, 10)),
                   "[0-9]+ terms were removed from gmt because they did not make the geneset.filter")
    expect_error(ActivePathways(dat, gmt, geneset.filter=c(0, NA)), NA)
    expect_error(ActivePathways(dat, gmt, geneset.filter=NULL), NA)
}) 
