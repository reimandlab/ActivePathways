# Prepare testing data
gmt <- read.GMT('test.gmt')
gmt_reac <- read.GMT('hsapiens_REAC_subset.gmt')
dat <- as.matrix(read.table('test_data.txt', header=TRUE, row.names='Gene'))
dat[is.na(dat)] <- 1
background <- makeBackground(gmt)

# filenames for cytoscape
CStag = "CS_files"

# Prepare testing data for scores_direction and constraints_vector
df <- read.table('test_data_rna_protein.tsv', header = TRUE, row.names = "gene", sep = '\t')
scores_test <- data.frame(row.names = rownames(df), rna = df$rna_pval, protein = df$protein_pval)
scores_test <- as.matrix(scores_test)
scores_test[is.na(scores_test)] <- 1
direction_test <- data.frame(row.names = rownames(df), rna = df$rna_log2fc, protein = df$protein_log2fc)
direction_test <- as.matrix(direction_test)
direction_test[is.na(direction_test)] <- 0
constraints_vector_test <- c(1,1)

# Run ActivePathways quickly
run_ap_short <- function(dat) ActivePathways(dat[,1, drop = F], gmt[1:3], cutoff=1, significant=1)
run_ap_short_contribution <- function(dat) ActivePathways(dat, gmt[1:3], cutoff=1, significant=1)
run_ap <- function(scores_test,direction_test,constraints_vector_test) ActivePathways(scores = scores_test,
                                                                                      gmt_reac, cutoff=1, significant=1,
                                                                                      scores_direction = direction_test,
                                                                                      constraints_vector = constraints_vector_test)



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
