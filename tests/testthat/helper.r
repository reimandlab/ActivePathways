# Prepare testing data
gmt <- read.GMT('test.gmt')
dat <- as.matrix(read.table('test_data.txt', header=TRUE, row.names='Gene'))
dat[is.na(dat)] <- 1
background <- makeBackground(gmt)

# filenames for cytoscape
CStag = "CS_files"

# Run ActivePathways quickly
run_ap_short <- function(dat) ActivePathways(dat[,1, drop = F], gmt[1:3], cutoff=1, significant=1)
run_ap_short_contribution <- function(dat) ActivePathways(dat, gmt[1:3], cutoff=1, significant=1)

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
