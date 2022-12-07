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

test_that("scores_direction and expected_direction have valid input",{
  
  dir_test <- direction_test
  dir_test[1,1] <- NA
  expect_error(run_ap(scores_test,dir_test,expected_direction_test), 'scores_direction may not contain missing values')
  
  
  dir_test <- direction_test
  dir_test[1,1] <- 'a'
  expect_error(run_ap(scores_test,dir_test,expected_direction_test), 'scores_direction must be numeric')
  
  
  dir_test <- direction_test
  rownames(dir_test) <- 1:length(direction_test[,1])
  expect_error(run_ap(scores_test,dir_test,expected_direction_test), 'scores_direction gene names must match scores genes')

  expected_dir <- c('a','b')
  expect_error(run_ap(scores_test,direction_test,expected_dir), 'expected_direction must be a numeric vector')
  
  expected_dir <- c(1,1,-1)
  expect_error(run_ap(scores_test,direction_test,expected_dir), 
               'expected_direction should have the same number of entries as columns in scores_direction')
  
  expected_dir <- c(1,-1)
  names(expected_dir) <- c("protein","rna")
  expect_error(run_ap(scores_test,direction_test,expected_dir), 
               'the expected_direction entries should match the order of scores_direction columns')
  
  
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

test_that("geneset_filter is a numeric vector of length 2", {
    expect_error(ActivePathways(dat, gmt, geneset_filter=1), 
                 "geneset_filter must be length 2")
    expect_error(ActivePathways(dat, gmt, geneset_filter=list(1,2)),
                 "geneset_filter must be a numeric vector")
    expect_error(ActivePathways(dat, gmt, geneset_filter=c('q', 2)),
                 "geneset_filter must be a numeric vector")
    expect_error(ActivePathways(dat, gmt, geneset_filter=c(1, -2)),
                 "geneset_filter limits must be positive")
    expect_error(ActivePathways(dat, gmt, geneset_filter=c(0, 0)), 
                 "No pathways in gmt made the geneset_filter", fixed=TRUE)
    expect_message(ActivePathways(dat, gmt, geneset_filter=c(NA, 10)),
                   "[0-9]+ terms were removed from gmt because they did not make the geneset_filter")
    expect_error(ActivePathways(dat, gmt, geneset_filter=c(0, NA)), NA)
    expect_error(ActivePathways(dat, gmt, geneset_filter=NULL), NA)
}) 

test_that("custom colors is a character vector that is equal in length to the number of columns in scores",{
  expect_error(ActivePathways(scores = dat, gmt = gmt, custom_colors = list("red","blue", "green")),
               "colors must be provided as a character vector",fixed = TRUE)  
  expect_error(ActivePathways(scores = dat, gmt = gmt, custom_colors = c("red","blue")),
               "incorrect number of colors is provided",fixed = TRUE)
  
  incorrect_color_names <- c("red","blue", "green")
  names(incorrect_color_names) <- c("promoter","lds",	"enhancer")
  expect_error(ActivePathways(scores = dat, gmt = gmt, custom_colors = incorrect_color_names),
               "names() of the custom colors vector should match the scores column names",fixed = TRUE) 
})

test_that("color palette is from the RColorBrewer package",{
  expect_error(ActivePathways(scores = dat, gmt = gmt, color_palette = "flamingo"),
               "palette must be from the RColorBrewer package",fixed = TRUE)
})

test_that("color palette and custom colors parameters are never specified together",{
  expect_error(ActivePathways(scores = dat, gmt = gmt, color_palette = "Pastel1", custom_colors = c("red","blue", "green")),
               "Both custom_colors and color_palette are provided. Specify only one of these parameters for node coloring.",fixed = TRUE)
})

test_that("color_integrated_only is a character vector of length 1",{
  expect_error(ActivePathways(scores = dat, gmt = gmt, color_integrated_only = list(1,2,3)),
               "color must be provided as a character vector",fixed = TRUE)
  expect_error(ActivePathways(scores = dat, gmt = gmt, color_integrated_only = c("red","blue")),
               "only a single color must be specified",fixed = TRUE)
})
