#' Merge a list or matrix of p-values
#'
#' @param scores Either a list of p-values or a matrix where each column is a test
#' @param method Method to merge p-valeus. See metap documentation for more info
#'
#' @return If scores is a list, returns a number. If scores is a matrix, returns
#'   a named list of p-values merged by row
merge_p_values <- function(scores, method=c("Fisher", "Brown", "logitp",
                                            "meanp", "sump", "sumz", "sumlog")) {
    if (ncol(scores) == 1) return (scores[, 1])

    method <- match.arg(method)
    if(method == "Fisher") method <- "sumlog"

    if (method == "Brown") {
        if (is.list(scores)) {
            stop("Brown's method cannot be used with a single list of p-values")
        }
        cov.matrix <- calculateCovariances(t(scores))
        return(apply(scores, 1, brownsMethod, cov.matrix=cov.matrix))
    }

    # Some metap function don't like p-values that are 0 or 1 so make them (0,1) to avoid errors
    scores <- apply(scores, c(1,2), function(x) if (x == 0) 0.0000001 else if (x==1) 0.9999999 else x)

    func <- function(x) get(method)(x)$p
    if (is.list(scores)) return(func(scores))
    return (apply(scores, 1, func))
}


#' Merge p-values using Brown's method
#'
#' Based on the R package EmpiricalBrownsMethod
#' https://github.com/IlyaLab/CombiningDependentPvaluesUsingEBM/blob/master/R/EmpiricalBrownsMethod/R/ebm.R
#' Only significant differences are the removal of extra_info and allowing a
#' pre-calculated covariance matrix
#'
#' @param p.values A vector of m p-values
#' @param data.matrix An m x n matrix representing m tests and n samples
#' @param cov.matrix A pre-calculated covariance matrix of data.matrix. More
#'   efficient when making many calls with the same data.matrix.
#'   Only one of data.matrix and cov.matrix must be given. If both are supplied,
#'   data.matrix is ignored
#' @return a p-value
brownsMethod <- function(p.values, data.matrix=NULL, cov.matrix=NULL) {
    if (missing(data.matrix) && missing(cov.matrix)) {
        stop ("Either data.matrix or cov.matrix must be supplied")
    }
    if (!(missing(data.matrix) || missing(cov.matrix))) {
        warning("Both data.matrix and cov.matrix were supplied. Ignoring data.matrix")
    }
    if (missing(cov.matrix)) cov.matrix <- calculateCovariances(data.matrix)

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
