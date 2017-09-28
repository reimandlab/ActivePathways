# Merge a list or matrix of p-values
# INPUT:
# scores:   Either a vector of p-values, or a numerical matrix where each column is a test   
# method:   method to merge p-values. See metap documentation for more infor on methods
# RETURNS:
# If scores is a list, returns a number
# If scores is a matrix, returns a named list of p-values merged by row
merge_p_values <- function(scores, method=c("Fisher", "Brown", "logtip", "meanp", "sump", "sumz", "sumlog")) {
    
    suppressPackageStartupMessages(metap)
    method <- match.arg(method)

    if (is.list(scores)) {
        if (method == "Brown") stop("Brown's method cannot be used with a single list of p-values")
        return (match.fun(method)(scores)$p)
    }

    if (method == "Brown") {
        corr.matrix <- calculateCovariances(t(scores))
        apply(scores, 1, empiricalBrownsMethod, corr.matrix=corr.matrix)
    }

    apply(scores, 1, match.fun(method)(scores)$p)
}



### BROWN'S METHOD ###
### Based on R package EmpiricalBrownsMethod
### https://github.com/IlyaLab/CombiningDependentPvaluesUsingEBM/blob/master/R/EmpiricalBrownsMethod/R/ebm.R
### The only significant differences are the removal of extra_info and allowing a pre-calculated correlation matrix

# INPUT
# data_matrix:      m x n matrix representing m tests and n samples
# p_values:         vector length m of p values
# cov.matrix:       pre-calculated covariance matrix of data_matrix. Saves time by
#                   not having to calculate the corvariance matrix for each call
brownsMethod <- function(data.matrix, p.values, cov.matrix=NULL) {
    #TODO: figure out how to handle NAs in p_values
    data.matrix <- na.exclude(data.matrix)

    if (is.null(cov.matrix)) cov.matrix <- calculateCovariances(data.matrix)

    N <- ncol(cov.matrix)
    expected <- 2 * N
    cov.sum <- 2 * sum(cov.matrix[lower.tri(cov.matrix, diag=FALSE)])
    var <- (4 * N) + cov.sum
    sf <- var / (2 * expected)

    df <- (2 * expected^2) / var
    if (df > 2 * N) {
        df <- 2 * N
        sf <- 1
    }
    
    x <- 2 * sum(-log(p.values), na.rm=TRUE)
    p.brown <- pchisq(df=df, q=x/sf, lower.tail=FALSE)
    p.brown
}

transformData <- function(dat) {
    dvm <- mean(dat, na.rm=TRUE)
    dvsd <- pop.sd(dat)
    s <- (dat - dvm) / dvsd
    distr <- ecdf(s)
    sapply(s, function(a) -2 * log(distr(a)))
}


calculateCovariances <- function(data.matrix) {
    transformed.data.matrix <- apply(data.matrix, 1, transformData)
    cov(transformed.data.matrix)
}

pop.var <- function(x) var(x, na.rm=TRUE) * (length(x) - 1) / length(x)
pop.sd <- function(x) sqrt(pop.var(x))
