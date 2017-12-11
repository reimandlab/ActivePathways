# Prepare testing data
gmt <- read.GMT('test.gmt')
dat <- as.matrix(read.table('test_data.txt', header=TRUE, row.names='Gene'))
dat[is.na(dat)] <- 1
background <- makeBackground(gmt)

# filenames for cytoscape
file.names <- paste(c('terms', 'groups', 'smallgmt'), '.txt', sep="")

# Run activeDriverPW quickly
run_adPW_short <- function(dat) activeDriverPW(dat, gmt[1:3], cutoff=1, significant=1, contribution=FALSE)
run_adPW_short_contribution <- function(dat) activeDriverPW(dat, gmt[1:3], cutoff=1, significant=1, contribution=TRUE)

# Data for testing enrichmentAnalysis
ea.gmt <- gmt[1:4]
ea.gmt[[1]]$genes <- c('PHC2', 'XPC', 'BLM')
ea.gmt[[2]]$genes <- c('HERC2', 'SP100', 'BLM')
ea.gmt[[3]]$genes <- c('HERC2', 'XPC')
ea.gmt[[4]]$genes <- c('XPC')
ea.genelist <- c('HERC2', 'SP100', 'BLM')
ea.genelist2 <- c('HERC2', letters, 'XPC')
ea.background <- makeBackground(ea.gmt)
ea.background2 <- c(ea.background, letters)

# Expectation to test if two lists contain the same items, ignoring order
expect_setequal <- function(actual, expected) {
    # Test that the sets from two objects are the same
    differences <- setdiff(actual,expected)
    sets_equal  <- length(differences) == 0
    message     <- paste("Sets not equal. First difference was:", differences[[1]])
    expect(sets_equal, message)
    invisible(actual)
}
