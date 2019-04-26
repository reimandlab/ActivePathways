#' Merge a list or matrix of p-values
#'
#' @param scores Either a list of p-values or a matrix where each column is a test
#' @param method Method to merge p-values. See 'methods' section below.
#'
#' @return If scores is a vector or list, returns a number. If scores is a
#'   matrix, returns a named list of p-values merged by row
#'
#' @section Methods:
#' There are multiple methods available that can be used to merge a list of
#'   p-values. The main methods encouraged by this package are:
#' \describe{
#'  \item{Fisher or sumlog}{Fisher's method assumes p-values are uniformly
#'  distributed and performs a chi-squared test on the statistic sum(-2 log(p)).
#'  This method is most appropriate when the columns in \code{scores} are
#'  independent.}
#'  \item{Brown}{Brown's method extends Fisher's method by accounting for the
#'  covariance in the columns of \code{scores}. It is more appropriate when the
#'  tests of significance used to create the columns in \code{scores} are not
#'  necessarily independent. Note that the "Brown" method cannot be used with a 
#'  single list of p-values. However, in this case Brown's method is identical 
#'  to Fisher's method.}
#' }
#' Other methods are also available. See \code{\link[metap]{metap-package}}
#' for more details
#'
#' @examples
#' merge_p_values(c(0.05, 0.09, 0.01))
#' merge_p_values(list(a=0.01, b=1, c=0.0015, d=0.025), method='meanp')
#' merge_p_values(matrix(data=c(0.03, 0.061, 0.48, 0.052), nrow=2), method='Brown')
#' 
#' @export
merge_p_values <- function(scores, method=c("Fisher", "Brown", "logitp",
                                            "meanp", "sump", "sumz", "sumlog")) {
    # Validation on scores
    if (is.list(scores)) scores <- unlist(scores, recursive=FALSE)
    if (!(is.vector(scores) || is.matrix(scores))) stop("scores must be a matrix or list")
    if (any(is.na(scores))) stop("scores may not contain missing values")
    if (!is.numeric(scores)) stop("scores must be numeric")
    if (any(scores < 0 | scores > 1)) stop("All values in scores must be in [0,1]")
        
    method <- match.arg(method)
    if(method == "Fisher") method <- "sumlog"

    if (is.vector(scores)) {
        if (method == "Brown") stop("Brown's method cannot be used with a single list of p-values")
        
        # Some metap functions don't like p-values that are 0 or 1 so make them (0, 1) to avoid errors
        scores <- sapply(scores, function(x) if (x == 0) 1e-16 else if (x==1) 1-1e-16 else x)
        func <- function(x) utils::getFromNamespace(method, 'metap')(x)$p
        return(func(scores))
    }
    
    # scores is a matrix
    if (ncol(scores) == 1) return (scores[, 1, drop=TRUE])
    
    if (method == "Brown") {
        cov.matrix <- calculateCovariances(t(scores))
        return(apply(scores, 1, brownsMethod, cov.matrix=cov.matrix))
    }
    
    scores <- apply(scores, c(1,2), function(x) if (x == 0) 1e-16 else if (x==1) 1-1e-16 else x)
    func <- function(x) utils::getFromNamespace(method, 'metap')(x)$p
    return (apply(scores, 1, func))
}


#' Merge p-values using Brown's method
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
        message("Both data.matrix and cov.matrix were supplied. Ignoring data.matrix")
    }
    if (missing(cov.matrix)) cov.matrix <- calculateCovariances(data.matrix)

    N <- ncol(cov.matrix)
    expected <- 2 * N
    cov.sum <- 2 * sum(cov.matrix[lower.tri(cov.matrix, diag=FALSE)])
    variance <- (4 * N) + cov.sum
    sf <- variance / (2 * expected)

    df <- (2 * expected^2) / variance
    if (df > 2 * N) {
        df <- 2 * N
        sf <- 1
    }

    x <- 2 * sum(-log(p.values), na.rm=TRUE)
    p.brown <- stats::pchisq(df=df, q=x/sf, lower.tail=FALSE)
    p.brown
}

transformData <- function(dat) {
    # If all values in dat are the same (equal to y), return dat. The covariance
    # matrix will be the zero matrix, and brown's method gives the p-value as y
    # Otherwise (dat - dmv) / dvsd is NaN and ecdf throws and error
    if (isTRUE(all.equal(min(dat), max(dat)))) return(dat)

    dvm <- mean(dat, na.rm=TRUE)
    dvsd <- pop.sd(dat)
    s <- (dat - dvm) / dvsd
    distr <- stats::ecdf(s)
    sapply(s, function(a) -2 * log(distr(a)))
}


calculateCovariances <- function(data.matrix) {
    transformed.data.matrix <- apply(data.matrix, 1, transformData)
    stats::cov(transformed.data.matrix)
}

pop.var <- function(x) stats::var(x, na.rm=TRUE) * (length(x) - 1) / length(x)
pop.sd <- function(x) sqrt(pop.var(x))
