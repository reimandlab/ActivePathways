#' Merge a list or matrix of p-values
#'
#' @param scores Either a list/vector of p-values or a matrix where each column is a test.
#' @param method Method to merge p-values. See 'methods' section below.
#' @param scores_direction Either a vector of log2 transformed fold-change values or a matrix where each column is a test. 
#' Must contain the same dimensions as the scores parameter. Datasets without directional information should be set to 0. 
#' @param constraints_vector  A numerical vector of +1 or -1 values corresponding to the user-defined
#'   directional relationship between the columns in scores_direction. Datasets without directional information should
#'   be set to 0. 
#'
#' @return If \code{scores} is a vector or list, returns a number. If \code{scores} is a
#'   matrix, returns a named list of p-values merged by row.
#'
#' @section Methods:
#' Four methods are available to merge a list of p-values:
#' \describe{
#'  \item{Fisher}{Fisher's method (default) assumes that p-values are uniformly
#'  distributed and performs a chi-squared test on the statistic sum(-2 log(p)).
#'  This method is most appropriate when the columns in \code{scores} are
#'  independent.}
#'  \item{Brown}{Brown's method extends Fisher's method by accounting for the
#'  covariance in the columns of \code{scores}. It is more appropriate when the
#'  tests of significance used to create the columns in \code{scores} are not
#'  necessarily independent. Note that the "Brown" method cannot be used with a 
#'  single list of p-values. However, in this case Brown's method is identical 
#'  to Fisher's method and should be used instead.}
#'  \item{Stouffer}{Stouffer's method assumes p-values are uniformly distributed
#'  and transforms p-values into a Z-score using the cumulative distribution function of a
#'  standard normal distribution. This method is appropriate when the columns in \code{scores}
#'   are independent.}
#'  \item{Strube}{Strube's method extends Stouffer's method by accounting for the 
#'  covariance in the columns of \code{scores}.}
#'   
#' }
#'
#' @examples
#'   merge_p_values(c(0.05, 0.09, 0.01))
#'   merge_p_values(list(a=0.01, b=1, c=0.0015, d=0.025), method='Fisher')
#'   merge_p_values(matrix(data=c(0.03, 0.061, 0.48, 0.052), nrow = 2), method='Brown')
#' 
#' @export
merge_p_values <- function(scores, method = "Fisher", scores_direction = NULL, 
                           constraints_vector = NULL) {
    
    ##### Validation #####
    # scores
    if (is.list(scores)) scores <- unlist(scores, recursive=FALSE)
    if (!(is.vector(scores) || is.matrix(scores))) stop("scores must be a matrix or vector")
    if (any(is.na(scores))) stop("scores cannot contain missing values, we recommend replacing NA with 1 or removing")
    if (!is.numeric(scores)) stop("scores must be numeric")
    if (any(scores < 0 | scores > 1)) stop("All values in scores must be in [0,1]")
    
    # method
    if (!method %in% c("Fisher", "Brown", "Stouffer","Strube")){
        stop("Only Fisher's, Brown's, Stouffer's and Strube's methods are currently supported")
    }
    
    # scores_direction and constraints_vector
    if (xor(!is.null(scores_direction),!is.null(constraints_vector))) stop("Both scores_direction and constraints_vector must be provided")
    if (!is.null(scores_direction) && !is.null(constraints_vector)){
        if (!(is.numeric(constraints_vector) && is.vector(constraints_vector))) stop("constraints_vector must be a numeric vector")
        if (any(!constraints_vector %in% c(1,-1,0))) stop("constraints_vector must contain the values: 1, -1 or 0")
        if (!(is.vector(scores_direction) || is.matrix(scores_direction))) stop("scores_direction must be a matrix or vector")
        if (!all(class(scores_direction) == class(scores))) stop("scores and scores_direction must be the same data type")
        if (any(is.na(scores_direction))) stop("scores_direction cannot contain missing values, we recommend replacing NA with 0 or removing")
        if (!is.numeric(scores_direction)) stop("scores_direction must be numeric")
        
        if (is.matrix(scores_direction)){
            if (any(!rownames(scores_direction) %in% rownames(scores))) stop ("scores_direction gene names must match scores genes")
            if (is.null(colnames(scores)) || is.null(colnames(scores_direction))) stop("column names must be provided to scores and scores_direction")
            if (any(!colnames(scores_direction) %in% colnames(scores))) stop("scores_direction column names must match scores column names")
            if (length(constraints_vector) != length(colnames(scores_direction))) stop("constraints_vector should have the same number of entries as columns in scores_direction")
            if (any(constraints_vector %in% 0) &&  !all(scores_direction[,constraints_vector %in% 0] == 0)) 
                stop("scores_direction entries must be set to 0's for columns that do not contain directional information")
            if (!is.null(names(constraints_vector))){
                if (!all.equal(names(constraints_vector), colnames(scores_direction), colnames(scores)) == TRUE){
                    stop("the constraints_vector entries should match the order of scores and scores_direction columns")
                }}}
        
        if (is.vector(scores_direction)){
            if (length(constraints_vector) != length(scores_direction)) stop("constraints_vector should have the same number of entries as scores_direction")
            if (length(scores_direction) != length(scores)) stop("scores_direction should have the same number of entries as scores")
            if (any(constraints_vector %in% 0) &&  !all(scores_direction[constraints_vector %in% 0] == 0)) 
                stop("scores_direction entries that do not contain directional information must be set to 0's")
            if (!is.null(names(constraints_vector))){
                if (!all.equal(names(constraints_vector), names(scores_direction), names(scores)) == TRUE){
                    stop("the constraints_vector entries should match the order of scores and scores_direction")
                }}}}
    
    ##### Merge P-values #####
    
    # Methods to merge p-values from a scores vector
    if (is.vector(scores)){
        if (method == "Brown" || method == "Strube") {
            stop("Brown's or Strube's method cannot be used with a single list of p-values")
        }
        
        # Convert zeroes to smallest available doubles
        scores <- sapply(scores, function(x) ifelse (x == 0, 1e-300, x))
        if (method == "Fisher"){
            p_fisher <- stats::pchisq(fishersMethod(scores, scores_direction,constraints_vector),
                                      2*length(scores), lower.tail = FALSE)
            return(p_fisher)
        } 
        if (method == "Stouffer"){
            p_stouffer <- 2*stats::pnorm(-1*abs(stouffersMethod(scores,scores_direction,constraints_vector)))
            return(p_stouffer)
        } 
    }
    
    # If scores is a matrix with one column, then no p-value merging can be done 
    if (ncol(scores) == 1) return (scores[, 1, drop=TRUE])
    
    # If scores is a matrix with multiple columns, apply the following methods
    scores <- apply(scores, c(1,2), function(x) ifelse (x == 0, 1e-300, x))
    
    if (method == "Fisher"){
        fisher_merged <- c()
        for(i in 1:length(scores[,1])) {
            p_fisher <- stats::pchisq(fishersMethod(scores[i,], scores_direction[i,],constraints_vector),
                                      2*length(scores[i,]), lower.tail = FALSE)
            fisher_merged <- c(fisher_merged,p_fisher)
        }
        names(fisher_merged) <- rownames(scores)
        return(fisher_merged)
    }
    if (method == "Brown") {
        cov_matrix <- calculateCovariances(t(scores))
        brown_merged <- brownsMethod(scores, cov_matrix = cov_matrix, scores_direction = scores_direction,
                                     constraints_vector = constraints_vector)
        return(brown_merged)
    }
    if (method == "Stouffer"){
        stouffer_merged <- c()
        for(i in 1:length(scores[,1])){
            p_stouffer <- 2*stats::pnorm(-1*abs(stouffersMethod(scores[i,], scores_direction[i,],constraints_vector)))
            stouffer_merged <- c(stouffer_merged,p_stouffer)
        }
        names(stouffer_merged) <- rownames(scores)
        return(stouffer_merged)
    }
    if (method == "Strube"){
        strube_merged <- strubesMethod(scores,scores_direction,constraints_vector)
        return(strube_merged)
    }
}


fishersMethod <- function(p_values, scores_direction, constraints_vector) {
    if (!is.null(scores_direction)){
        directionality <- constraints_vector * scores_direction/abs(scores_direction)
        p_values_directional <- p_values[!is.na(directionality)]
        chisq_directional <- abs(-2 * sum(log(p_values_directional)*directionality[!is.na(directionality)]))
        chisq_nondirectional <- abs(-2 * sum(log(p_values[is.na(directionality)])))
        chisq_values <- c(chisq_directional, chisq_nondirectional)
    } else {
        chisq_values <- -2*log(p_values)
    }
    sum(chisq_values)
}


#' Merge p-values using the Brown's method. If scores_direction and constraints_vector parameters
#'  are provided this method is called DPM 
#'
#' @param p_values A matrix of m x n p-values.
#' @param data_matrix An m x n matrix representing m tests and n samples. NA's are not allowed.
#' @param cov_matrix A pre-calculated covariance matrix of \code{data_matrix}. This is more
#'   efficient when making many calls with the same data_matrix.
#'   Only one of \code{data_matrix} and \code{cov_matrix} must be given. If both are supplied,
#'   \code{data_matrix} is ignored.
#' @param scores_direction A matrix of log2 fold-change values. Datasets without directional information should be set to 0. 
#' @param constraints_vector  A numerical vector of +1 or -1 values corresponding to the user-defined
#'   directional relationship between columns in scores_direction. Datasets without directional information should
#'   be set to 0.
#' @return A p-value vector representing the merged significance of multiple p-values.
#' @export

# Based on the R package EmpiricalBrownsMethod
# https://github.com/IlyaLab/CombiningDependentPvaluesUsingEBM/blob/master/R/EmpiricalBrownsMethod/R/ebm.R
# Only significant differences are the removal of extra_info and allowing a
# pre-calculated covariance matrix
# 
brownsMethod <- function(p_values, data_matrix = NULL, cov_matrix = NULL, scores_direction,
                         constraints_vector) {
    if (missing(data_matrix) && missing(cov_matrix)) {
        stop ("Either data_matrix or cov_matrix must be supplied")
    }
    if (!(missing(data_matrix) || missing(cov_matrix))) {
        message("Both data_matrix and cov_matrix were supplied. Ignoring data_matrix")
    }
    if (missing(cov_matrix)) cov_matrix <- calculateCovariances(data_matrix)
    
    N <- ncol(cov_matrix)
    expected <- 2 * N
    cov_sum <- 2 * sum(cov_matrix[lower.tri(cov_matrix, diag=FALSE)])
    var <- (4 * N) + cov_sum
    sf <- var / (2 * expected)
    
    df <- (2 * expected^2) / var
    if (df > 2 * N) {
        df <- 2 * N
        sf <- 1
    }
    
    # Acquiring the unadjusted chi-squared value from Fisher's method
    fisher_chisq <- c()
    for(i in 1:length(p_values[,1])) {
        fisher_chisq <- c(fisher_chisq, fishersMethod(p_values[i,], scores_direction[i,],constraints_vector))
    }
    
    # Adjusted p-value
    p_brown <- stats::pchisq(df=df, q=fisher_chisq/sf, lower.tail=FALSE)
    names(p_brown) <- rownames(p_values)
    p_brown
}


stouffersMethod <- function (p_values, scores_direction,constraints_vector){
    k = length(p_values)
    if (!is.null(scores_direction)){
        directionality <- constraints_vector * scores_direction/abs(scores_direction)
        p_values_directional <- p_values[!is.na(directionality)]
        z_directional <- abs(sum(stats::qnorm(p_values_directional/2)*directionality[!is.na(directionality)]))
        z_nondirectional <- abs(sum(stats::qnorm(p_values[is.na(directionality)]/2)))
        z_values <- c(z_directional, z_nondirectional)
    } else {
        z_values <- stats::qnorm(p_values/2)
    }
    sum(z_values)/sqrt(k)
}


strubesMethod <- function (p_values, scores_direction, constraints_vector){
    # Acquiring the unadjusted z-value from Stouffer's method
    stouffer_z <- c()
    for(i in 1:length(p_values[,1])){
        stouffer_z <- c(stouffer_z,stouffersMethod(p_values[i,], scores_direction[i,],constraints_vector))
    }
    
    # Correlation matrix
    cor_mtx <- stats::cor(p_values, use = "complete.obs")
    cor_mtx[is.na(cor_mtx)] <- 0
    cor_mtx <- abs(cor_mtx)
    
    # Adjusted p-value
    k = length(p_values[1,])
    adjusted_z <- stouffer_z * sqrt(k) / sqrt(sum(cor_mtx))
    p_strube <- 2*stats::pnorm(-1*abs(adjusted_z))
    names(p_strube) <- rownames(p_values)
    p_strube
}


transformData <- function(dat) {
    # If all values in dat are the same (equal to y), return dat. The covariance
    # matrix will be the zero matrix, and brown's method gives the p-value as y
    # Otherwise (dat - dmv) / dvsd is NaN and ecdf throws an error
    if (isTRUE(all.equal(min(dat), max(dat)))) return(dat)
    
    dvm <- mean(dat, na.rm=TRUE)
    dvsd <- pop.sd(dat)
    s <- (dat - dvm) / dvsd
    distr <- stats::ecdf(s)
    sapply(s, function(a) -2 * log(distr(a)))
}


calculateCovariances <- function(data_matrix) {
    transformed_data_matrix <- apply(data_matrix, 1, transformData)
    stats::cov(transformed_data_matrix)
}

pop.var <- function(x) stats::var(x, na.rm=TRUE) * (length(x) - 1) / length(x)
pop.sd <- function(x) sqrt(pop.var(x))
