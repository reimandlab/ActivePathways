context("merge_p_values function")

test_list <- list(a=0.01, b=0.06, c=0.8, d=0.0001, e=0, f=1)
test_matrix <- matrix(c(0.01, 0.06, 0.08, 0.0001, 0, 1), ncol=2)

comparison_list <- test_list
comparison_list[[5]] = 1e-300

test_matrix = matrix(unlist(test_list), ncol = 2)
comparison_matrix = matrix(unlist(comparison_list), ncol = 2)

test_that("scores is a numeric matrix or list with valid p-values", {

    expect_error(merge_p_values(unlist(test_list), "Fisher"), NA)
    expect_error(merge_p_values(test_list, "Brown"), 
                 "Brown's or Strube's method cannot be used with a single list of p-values")
    expect_error(merge_p_values(unlist(test_list), "Brown"), 
                 "Brown's or Strube's method cannot be used with a single list of p-values")
    
    
    
    test_list[[1]] <- -0.1
    expect_error(merge_p_values(test_list), 'All values in scores must be in [0,1]', fixed=TRUE)
    test_list[[1]] <- 1.1
    expect_error(merge_p_values(test_list), 'All values in scores must be in [0,1]', fixed=TRUE)
    test_list[[1]] <- NA
    expect_error(merge_p_values(test_list), 'scores may not contain missing values')
    test_list[[1]] <- 'c'
    expect_error(merge_p_values(test_list), 'scores must be numeric')
    
    
    test_matrix[1, 1] <- NA
    expect_error(merge_p_values(test_matrix), 'scores may not contain missing values')
    test_matrix[1, 1] <- -0.1
    expect_error(merge_p_values(test_matrix), "All values in scores must be in [0,1]", fixed=TRUE)
    test_matrix[1, 1] <- 1.1
    expect_error(merge_p_values(test_matrix), "All values in scores must be in [0,1]", fixed=TRUE)
    test_matrix[1, 1] <- 'a'
    expect_error(merge_p_values(test_matrix), 'scores must be numeric')
    
    
})



test_that("Merged p-values are correct", {
	
	this_tolerance = 1e-7
	answer1 = c(1.481551e-05, 4.167534e-299, 9.785148e-01)
	answer2 = c(2.52747e-05, 0.00000e+00, 9.73873e-01)
	answer3 = 7.147579e-296
    
    expect_equal(merge_p_values(test_matrix, "Fisher"), answer1, tolerance = this_tolerance)

    expect_equal(merge_p_values(test_matrix, "Brown"), answer2, tolerance = this_tolerance)
    
    expect_equal(merge_p_values(test_matrix[, 1, drop=FALSE], "Fisher"), test_matrix[, 1, drop=TRUE])
    expect_equal(merge_p_values(test_matrix[, 1, drop=FALSE], "Brown"), test_matrix[, 1, drop=TRUE])

    expect_equal(merge_p_values(test_list, "Fisher"), answer3, tolerance = this_tolerance)
})


test_matrix <- matrix(c(0.01, 0.06, 0.08, 0.0001, 0, 1), ncol=2)
test_direction_matrix <- matrix(c(1,-1,1,-1,-1,1), ncol=2)
expected_dir <- c(1,1)

colnames(test_matrix) <- c("RNA", "Protein") 
rownames(test_matrix) <- c("TP53", "CHRNA1","PTEN")

colnames(test_direction_matrix) <- colnames(test_matrix)
rownames(test_direction_matrix) <- rownames(test_matrix)

test_that("scores_direction and expected_direction are valid", {
      
      expect_error(merge_p_values(test_matrix, "Fisher",test_direction_matrix),'Both scores_direction and expected_direction must be provided')
      expect_error(merge_p_values(test_matrix, "Fisher",expected_direction = expected_dir),'Both scores_direction and expected_direction must be provided')
      expect_error(merge_p_values(test_matrix, "Fisher", test_direction_matrix, c(1,"b")), 'expected_direction must be a numeric vector')
      expect_error(merge_p_values(test_matrix, "Fisher", test_direction_matrix, c(1,0)), "scores_direction entries must be set to 0's for columns that do not contain directional information")
      expect_error(merge_p_values(test_matrix, "Fisher", test_direction_matrix, c(1,5)), "expected_direction must contain the values: 1, -1 or 0")
      
      test_dir <- as.vector(test_direction_matrix)
      expect_error(merge_p_values(test_matrix, "Fisher", test_dir, c(1,1)), 'scores and scores_direction must be the same data type')
      
      test_dir <- test_direction_matrix
      test_dir[1,1] <- NA
      expect_error(merge_p_values(test_matrix, "Fisher", test_dir, c(1,1)), 'scores_direction may not contain missing values')
      
      test_dir <- test_direction_matrix
      test_dir[1,1] <- 'a'
      expect_error(merge_p_values(test_matrix, "Fisher", test_dir, c(1,1)), 'scores_direction must be numeric')
      
      test_dir <- test_direction_matrix
      colnames(test_dir) <- NULL
      expect_error(merge_p_values(test_matrix, "Fisher", test_dir, c(1,1)), 'column names must be provided to scores and scores_direction')
      
      test_m <- test_matrix
      rownames(test_m) <- c("TP53", "GENE2", "GENE3")
      expect_error(merge_p_values(test_m, "Fisher", test_direction_matrix, c(1,1)), 'scores_direction gene names must match scores genes')
      
      test_m <- test_matrix
      colnames(test_m) <- c("RNA","Mutation")
      expect_error(merge_p_values(test_m, "Fisher", test_direction_matrix, c(1,1)),
                   'scores_direction column names must match scores column names')
      
      expect_error(merge_p_values(test_matrix, "Fisher", test_direction_matrix, c(1,1,-1)),
                   'expected_direction should have the same number of entries as columns in scores_direction')
      
      names(expected_dir) <- c("Protein","RNA")
      expect_error(merge_p_values(test_matrix, "Fisher", test_direction_matrix, expected_dir),
                   'the expected_direction entries should match the order of scores and scores_direction columns')
})

