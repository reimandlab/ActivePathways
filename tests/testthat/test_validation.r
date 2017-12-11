context("Validation on parameters")


test_that("scores is a numeric matrix with valid p-values", {
  dat2 <- dat
  dat2[1, 1] <- 'a'
  expect_error(run_adPW_short(dat2), 'scores must be a numeric matrix')

  dat2 <- dat
  dat2[1, 1] <- NA
  expect_error(run_adPW_short(dat2), 'scores may not contain missing values')

  dat2[1, 1] <- -0.1
  expect_error(run_adPW_short(dat2), "All values in scores must be in [0,1]", fixed=TRUE)

  dat2[1, 1] <- 1.1
  expect_error(run_adPW_short(dat2), "All values in scores must be in [0,1]", fixed=TRUE)

  dat2[1, 1] <- 1
  expect_error(run_adPW_short(dat2), NA)

  dat2[1, 1] <- 0
  expect_error(run_adPW_short(dat2), NA)
})

test_that("significant is valid", {
    expect_error(activeDriverPW(dat, gmt, significant=-0.1),
                 "significant must be a value in [0,1]", fixed=TRUE)
    expect_error(activeDriverPW(dat, gmt, significant = 1.1),
                 "significant must be a value in [0,1]", fixed=TRUE)
    expect_warning(activeDriverPW(dat, gmt, significant=0, return.all=TRUE),
                   "No significant terms were found")
    expect_error(activeDriverPW(dat, gmt, significant=1, return.all=TRUE), NA)
})


test_that("cutoff is valid", {
    expect_error(activeDriverPW(dat, gmt, cutoff=-0.1),
                 "cutoff must be a value in [0,1]", fixed=TRUE)
    expect_error(activeDriverPW(dat, gmt, cutoff = 1.1),
                 "cutoff must be a value in [0,1]", fixed=TRUE)
    expect_error(activeDriverPW(dat, gmt, cutoff=0, return.all=TRUE),
                 "No genes made the cutoff", fixed=TRUE)
    expect_error(activeDriverPW(dat, gmt, cutoff=1, return.all=TRUE), NA)
})


test_that("background is a character vector", {
    error_msg <- "background must be a character vector"
    expect_error(activeDriverPW(dat, gmt, background=c(1,5,2)), error_msg)
    expect_error(activeDriverPW(dat, gmt, background=matrix(c('a', 'b', 'c', 'd'), 2)), error_msg)
})


test_that("genes not found in background are removed", {
    expect_message(activeDriverPW(dat, gmt, background=rownames(dat)[-(1:10)], significant=1, cutoff=1),
                   "10 rows were removed from scores because they are not found in the background")
    expect_error(activeDriverPW(dat, gmt, background='qwerty'),
                 "scores does not contain any genes in the background")
})
