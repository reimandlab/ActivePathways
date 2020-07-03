context("merge_p_values function")

test_that("scores is a numeric matrix or list with valid p-values", {
    test.list <- list(a=0.01, b=0.06, c=0.8, d=0.0001, e=0, f=1)

    expect_error(merge_p_values(unlist(test.list), "Fisher"), NA)
    expect_error(merge_p_values(test.list, "Brown"), 
                 "Brown's method cannot be used with a single list of p-values")
    expect_error(merge_p_values(unlist(test.list), "Brown"), 
                 "Brown's method cannot be used with a single list of p-values")
    
    test.list[[1]] <- -0.1
    expect_error(merge_p_values(test.list), 'All values in scores must be in [0,1]', fixed=TRUE)
    test.list[[1]] <- 1.1
    expect_error(merge_p_values(test.list), 'All values in scores must be in [0,1]', fixed=TRUE)
    test.list[[1]] <- NA
    expect_error(merge_p_values(test.list), 'scores may not contain missing values')
    test.list[[1]] <- 'c'
    expect_error(merge_p_values(test.list), 'scores must be numeric')
    
    
    test.matrix <- matrix(c(0.01, 0.06, 0.08, 0.0001, 0, 1), ncol=2)
    test.matrix[1, 1] <- NA
    expect_error(merge_p_values(test.matrix), 'scores may not contain missing values')
    test.matrix[1, 1] <- -0.1
    expect_error(merge_p_values(test.matrix), "All values in scores must be in [0,1]", fixed=TRUE)
    test.matrix[1, 1] <- 1.1
    expect_error(merge_p_values(test.matrix), "All values in scores must be in [0,1]", fixed=TRUE)
    test.matrix[1, 1] <- 'a'
    expect_error(merge_p_values(test.matrix), 'scores must be numeric')
})

test_that("Merged p-values are correct", {
    test.list <- list(a=0.01, b=0.06, c=0.8, d=0.0001, e=0, f=1)
    comparison.list <- c(a=0.01, b=0.06, c=0.8, d=0.0001, e=1e-300, f=1)
    
    test.matrix <- matrix(c(0.01, 0.06, 0.08, 0.0001, 0, 1), ncol=2)
    comparison.matrix <- matrix(c(0.01, 0.06, 0.08, 0.0001, 1e-300, 1), ncol=2)
    
    expect_equal(merge_p_values(test.matrix, "Fisher"), apply(comparison.matrix, 1, 
    		function(x) metap::sumlog(x)$p))

    expect_equal(merge_p_values(test.matrix, "Brown"), 
			apply(test.matrix, 1, function(x) 
            	EmpiricalBrownsMethod::empiricalBrownsMethod(t(test.matrix), x, FALSE)))
    
    expect_equal(merge_p_values(test.matrix[, 1, drop=FALSE], "Fisher"), test.matrix[, 1, drop=TRUE])
    expect_equal(merge_p_values(test.matrix[, 1, drop=FALSE], "Brown"), test.matrix[, 1, drop=TRUE])
    expect_equal(merge_p_values(test.list, "Fisher"), metap::sumlog(comparison.list)$p)
})