# Prepare testing data
gmt <- read.GMT('test.gmt')
gmt_hox <- read.GMT('hsapiens_REAC_subset.gmt')
dat <- as.matrix(read.table('test_data.txt', header=TRUE, row.names='Gene'))
dat[is.na(dat)] <- 1
background <- makeBackground(gmt)

# filenames for cytoscape
CStag = "CS_files"

# Prepare testing data for scores_direction and expected_direction
df_hox <- read.table('test_hoxa10as_data.tsv', header = TRUE, sep = '\t')
hox_score <- data.frame(row.names = df_hox[,1], overexpression = df_hox[,3], knockdown = df_hox[,5])
hox_score <- as.matrix(hox_score)
hox_score[is.na(hox_score)] <- 1
hox_direction <- data.frame(row.names = df_hox[,1], overexpression = df_hox[,2], knockdown = df_hox[,4])
hox_direction <- as.matrix(hox_direction)
hox_direction[is.na(hox_direction)] <- 1
hox_expected_direction <- c(-1,1)

# Run ActivePathways quickly
run_ap_short <- function(dat) ActivePathways(dat[,1, drop = F], gmt[1:3], cutoff=1, significant=1)
run_ap_short_contribution <- function(dat) ActivePathways(dat, gmt[1:3], cutoff=1, significant=1)
run_ap <- function(hox_score,hox_direction,hox_expected_direction) ActivePathways(scores = hox_score,
                                                                                  gmt_hox, cutoff=1, significant=1,
                                                                                  scores_direction = hox_direction,
                                                                                  expected_direction = hox_expected_direction)



# Data for testing enrichmentAnalysis
ea_gmt <- gmt[1:4]
ea_gmt[[1]]$genes <- c('PHC2', 'XPC', 'BLM')
ea_gmt[[2]]$genes <- c('HERC2', 'SP100', 'BLM')
ea_gmt[[3]]$genes <- c('HERC2', 'XPC')
ea_gmt[[4]]$genes <- c('XPC')
ea_genelist <- c('HERC2', 'SP100', 'BLM')
ea_genelist2 <- c('HERC2', letters, 'XPC')
ea_background <- makeBackground(ea_gmt)
ea_background2 <- c(ea_background, letters)

# Expectation to test if two lists contain the same items, ignoring order
expect_setequal <- function(actual, expected) {
    # Test that the sets from two objects are the same
    differences <- setdiff(actual,expected)
    sets_equal  <- length(differences) == 0
    message     <- paste("Sets not equal. First difference was:", differences[[1]])
    expect(sets_equal, message)
    invisible(actual)
}
