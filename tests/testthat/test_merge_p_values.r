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

