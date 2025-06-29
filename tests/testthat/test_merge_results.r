context("merge_results function")

test_matrix <- matrix(c(0.01, 0.06, 0.08, 0.0001, 0, 1), ncol=2)
test_direction_matrix <- matrix(c(1, -1, 1, -1, -1, 1), ncol=2)
constraints_vector <- c(1, 1)

colnames(test_matrix) <- c("rna", "protein") 
rownames(test_matrix) <- c("TP53", "CHRNA1", "PTEN")

colnames(test_direction_matrix) <- colnames(test_matrix)
rownames(test_direction_matrix) <- rownames(test_matrix)

enriched_pathways <- ActivePathways(scores=test_matrix, gmt_reac, cutoff=1, significant=1)
enriched_pathways_directional <- run_ap(test_matrix, test_direction_matrix, constraints_vector)

test_that("merge_results inputs are valid", {

    # Check that test parameter is not NULL or length 0
    expect_error(merge_results(enriched_pathways, enriched_pathways_directional, "hsapiens_REAC_subset.gmt", tests=c()),"Tests parameter must be provided (e.g., c('rna', 'protein', 'combined')) or NULL", fixed=TRUE)
    expect_error(merge_results(enriched_pathways, enriched_pathways_directional, "hsapiens_REAC_subset.gmt", tests=NULL),"Tests parameter must be provided (e.g., c('rna', 'protein', 'combined')) or NULL", fixed=TRUE)

    # Check if all tests (except 'combined') exist in both dataframes with "Genes_" prefix
    expect_error(merge_results(enriched_pathways, enriched_pathways_directional, "hsapiens_REAC_subset.gmt", tests=c("methylation", "protein", "combined")),"All tests (except 'combined') must exist as columns with 'Genes_' prefix in both enriched_pathways and enriched_pathways_directional", fixed=TRUE)
    expect_error(merge_results(enriched_pathways[,1:3], enriched_pathways_directional, "hsapiens_REAC_subset.gmt", tests=c("rna", "protein", "combined")),"All tests (except 'combined') must exist as columns with 'Genes_' prefix in both enriched_pathways and enriched_pathways_directional", fixed=TRUE)
    expect_error(merge_results(enriched_pathways, enriched_pathways_directional[,1:3], "hsapiens_REAC_subset.gmt", tests=c("rna", "protein", "combined")),"All tests (except 'combined') must exist as columns with 'Genes_' prefix in both enriched_pathways and enriched_pathways_directional", fixed=TRUE)

    # Check that col_colors is the same length as tests
    expect_error(merge_results(enriched_pathways, enriched_pathways_directional, "hsapiens_REAC_subset.gmt", tests=c("rna", "protein", "combined"), col_colors=c("#FF0000", "#00FF00")),"col_colors must have the same length as tests", fixed=TRUE)
    expect_error(merge_results(enriched_pathways, enriched_pathways_directional, "hsapiens_REAC_subset.gmt", tests=c("rna", "protein", "combined"), col_colors=c("green")),"col_colors must have the same length as tests", fixed=TRUE)
   
})
