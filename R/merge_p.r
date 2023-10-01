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
#' Eight methods are available to merge a list of p-values:
#' \describe{
#'  \item{Fisher}{Fisher's method (default) assumes that p-values are uniformly
#'  distributed and performs a chi-squared test on the statistic sum(-2 log(p)).
#'  This method is most appropriate when the columns in \code{scores} are
#'  independent.}
#'  \item{Fisher_directional}{Fisher's method modification that allows for 
#'  directional information to be incorporated with the \code{scores_direction}
#'  and \code{constraints_vector} parameters.}
#'  \item{Brown}{Brown's method extends Fisher's method by accounting for the
#'  covariance in the columns of \code{scores}. It is more appropriate when the
#'  tests of significance used to create the columns in \code{scores} are not
#'  necessarily independent. Note that the "Brown" method cannot be used with a 
#'  single list of p-values. However, in this case Brown's method is identical 
#'  to Fisher's method and should be used instead.}
#'  \item{DPM}{DPM extends Brown's method by incorporating directional information
#'  using the \code{scores_direction} and \code{constraints_vector} parameters.}
#'  \item{Stouffer}{Stouffer's method assumes p-values are uniformly distributed
#'  and transforms p-values into a Z-score using the cumulative distribution function of a
#'  standard normal distribution. This method is appropriate when the columns in \code{scores}
#'   are independent.}
#'  \item{Stouffer_directional}{Stouffer's method modification that allows for 
#'  directional information to be incorporated with the \code{scores_direction}
#'  and \code{constraints_vector} parameters.}
#'  \item{Strube}{Strube's method extends Stouffer's method by accounting for the 
#'  covariance in the columns of \code{scores}.}
#'  \item{Strube_directional}{Strube's method modification that allows for 
#'  directional information to be incorporated with the \code{scores_direction}
#'  and \code{constraints_vector} parameters.}
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
  
  # scores_direction and constraints_vector
  if (xor(!is.null(scores_direction),!is.null(constraints_vector))) stop("Both scores_direction and constraints_vector must be provided")
  if (!is.null(scores_direction) && !is.null(constraints_vector)){
    if (!(is.numeric(constraints_vector) && is.vector(constraints_vector))) stop("constraints_vector must be a numeric vector")
    if (any(!constraints_vector %in% c(1,-1,0))) stop("constraints_vector must contain the values: 1, -1 or 0")
    if (!(is.vector(scores_direction) || is.matrix(scores_direction))) stop("scores_direction must be a matrix or vector")
    if (!all(class(scores_direction) == class(scores))) stop("scores and scores_direction must be the same data type")
    if (any(is.na(scores_direction))) stop("scores_direction cannot contain missing values, we recommend replacing NA with 0 or removing")
    if (!is.numeric(scores_direction)) stop("scores_direction must be numeric")
    if (method %in% c("Fisher","Brown","Stouffer","Strube")) stop("Only DPM, Fisher_directional, Stouffer_directional, and Strube_directional methods support directional integration")
    
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
  
  # method
  if (!method %in% c("Fisher", "Fisher_directional", "Brown", "DPM", "Stouffer", "Stouffer_directional", "Strube", "Strube_directional")){
    stop("Only Fisher, Brown, Stouffer and Strube methods are currently supported for non-directional analysis. 
             And only DPM, Fisher_directional, Stouffer_directional, and Strube_directional are supported for directional analysis")
  }
  if (method %in% c("Fisher_directional", "DPM", "Stouffer_directional", "Strube_directional") & 
      is.null(scores_direction)){
    stop("scores_direction and constraints_vector must be provided for directional analyses")
  }
  
  
  ##### Merge P-values #####
  
  # Methods to merge p-values from a scores vector
  if (is.vector(scores)){
    if (method == "Brown" || method == "Strube" || method == "DPM" || method == "Strube_directional") {
      stop("Brown's, DPM, Strube's, and Strube_directional methods cannot be used with a single list of p-values")
    }
    
    # Convert 0 or very small p-values to 1e-300
    if(min(scores) < 1e-300){  
      message(paste('warning: p-values smaller than ', 1e-300, ' are replaced with ', 1e-300))
      scores <- sapply(scores, function(x) ifelse (x < 1e-300, 1e-300, x))
    }
    
    if (method == "Fisher"){
      p_fisher <- stats::pchisq(fishersMethod(scores),2*length(scores), lower.tail = FALSE)
      return(p_fisher)
    }
    
    if (method == "Fisher_directional"){
      p_fisher <- stats::pchisq(fishersDirectional(scores, scores_direction,constraints_vector),
                                2*length(scores), lower.tail = FALSE)
      return(p_fisher)
    } 
    
    if (method == "Stouffer"){
      p_stouffer <- 2*stats::pnorm(-1*abs(stouffersMethod(scores)))
      return(p_stouffer)
    } 
    if (method == "Stouffer_directional"){
      p_stouffer <- 2*stats::pnorm(-1*abs(stouffersDirectional(scores,scores_direction,constraints_vector)))
      return(p_stouffer)
    } 
  }
  
  # If scores is a matrix with one column, then no p-value merging can be done 
  if (ncol(scores) == 1) return (scores[, 1, drop=TRUE])
  
  # If scores is a matrix with multiple columns, apply the following methods
  if(min(scores) < 1e-300){
    message(paste('warning: p-values smaller than ', 1e-300, ' are replaced with ', 1e-300))
    scores <- apply(scores, c(1,2), function(x) ifelse (x < 1e-300, 1e-300, x))
  }
  
  if (method == "Fisher"){
    fisher_merged <- c()
    for(i in 1:length(scores[,1])) {
      p_fisher <- stats::pchisq(fishersMethod(scores[i,]), 2*length(scores[i,]), lower.tail = FALSE)
      fisher_merged <- c(fisher_merged, p_fisher)
    }
    names(fisher_merged) <- rownames(scores)
    return(fisher_merged)
  }
  if (method == "Fisher_directional"){
    fisher_merged <- c()
    for(i in 1:length(scores[,1])) {
      p_fisher <- stats::pchisq(fishersDirectional(scores[i,], scores_direction[i,], constraints_vector),
                                2*length(scores[i,]), lower.tail = FALSE)
      fisher_merged <- c(fisher_merged,p_fisher)
    }
    names(fisher_merged) <- rownames(scores)
    return(fisher_merged)
  }
  if (method == "Brown") {
    cov_matrix <- calculateCovariances(t(scores))
    brown_merged <- brownsMethod(scores, cov_matrix = cov_matrix)
    return(brown_merged)
  }
  if (method == "DPM") {
    cov_matrix <- calculateCovariances(t(scores))
    dpm_merged <- DPM(scores, cov_matrix = cov_matrix, scores_direction = scores_direction,
                      constraints_vector = constraints_vector)
    return(dpm_merged)
  }
  if (method == "Stouffer"){
    stouffer_merged <- c()
    for(i in 1:length(scores[,1])){
      p_stouffer <- 2*stats::pnorm(-1*abs(stouffersMethod(scores[i,])))
      stouffer_merged <- c(stouffer_merged,p_stouffer)
    }
    names(stouffer_merged) <- rownames(scores)
    return(stouffer_merged)
  }
  if (method == "Stouffer_directional"){
    stouffer_merged <- c()
    for(i in 1:length(scores[,1])){
      p_stouffer <- 2*stats::pnorm(-1*abs(stouffersDirectional(scores[i,], scores_direction[i,],constraints_vector)))
      stouffer_merged <- c(stouffer_merged,p_stouffer)
    }
    names(stouffer_merged) <- rownames(scores)
    return(stouffer_merged)
  }
  if (method == "Strube"){
    strube_merged <- strubesMethod(scores)
    return(strube_merged)
  }
  if (method == "Strube_directional"){
    strube_merged <- strubesDirectional(scores,scores_direction,constraints_vector)
    return(strube_merged)
  }
}


fishersMethod <- function(p_values) {
  chisq_values <- -2*log(p_values)
  sum(chisq_values)
}


fishersDirectional <- function(p_values, scores_direction, constraints_vector) {
  # Sum the directional chi-squared values
  directionality <- constraints_vector * scores_direction/abs(scores_direction)
  p_values_directional <- p_values[!is.na(directionality)]
  chisq_directional <- abs(-2 * sum(log(p_values_directional)*directionality[!is.na(directionality)]))
  
  # Sum the non-directional chi-squared values
  chisq_nondirectional <- abs(-2 * sum(log(p_values[is.na(directionality)])))
  
  # Combine both
  sum(c(chisq_directional, chisq_nondirectional))
}


#' Merge p-values using the Brown's method.
#'
#' @param p_values A matrix of m x n p-values.
#' @param data_matrix An m x n matrix representing m tests and n samples. NA's are not allowed.
#' @param cov_matrix A pre-calculated covariance matrix of \code{data_matrix}. This is more
#'   efficient when making many calls with the same data_matrix.
#'   Only one of \code{data_matrix} and \code{cov_matrix} must be given. If both are supplied,
#'   \code{data_matrix} is ignored.
#' @return A p-value vector representing the merged significance of multiple p-values.
#' @export

# Based on the R package EmpiricalBrownsMethod
# https://github.com/IlyaLab/CombiningDependentPvaluesUsingEBM/blob/master/R/EmpiricalBrownsMethod/R/ebm.R
# Only significant differences are the removal of extra_info and allowing a
# pre-calculated covariance matrix
# 
brownsMethod <- function(p_values, data_matrix = NULL, cov_matrix = NULL) {
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
  
  # Acquiring the unadjusted chi-squared values from Fisher's method
  fisher_chisq <- c()
  for(i in 1:length(p_values[,1])) {
    fisher_chisq <- c(fisher_chisq, fishersMethod(p_values[i,]))
  }
  
  # Adjusted p-value
  p_brown <- stats::pchisq(df=df, q=fisher_chisq/sf, lower.tail=FALSE)
  names(p_brown) <- rownames(p_values)
  p_brown
}


#' Merge p-values using the DPM method.
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


DPM <- function(p_values, data_matrix = NULL, cov_matrix = NULL,
                scores_direction, constraints_vector) {
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
    fisher_chisq <- c(fisher_chisq, fishersDirectional(p_values[i,], scores_direction[i,],constraints_vector))
  }
  
  # Adjusted p-value
  p_dpm <- stats::pchisq(df=df, q=fisher_chisq/sf, lower.tail=FALSE)
  names(p_dpm) <- rownames(p_values)
  p_dpm
}


stouffersMethod <- function (p_values){
  k = length(p_values)
  z_values <- stats::qnorm(p_values/2)
  sum(z_values)/sqrt(k)
}


stouffersDirectional <- function (p_values, scores_direction, constraints_vector){
  k = length(p_values)
  
  # Sum the directional z-values
  directionality <- constraints_vector * scores_direction/abs(scores_direction)
  p_values_directional <- p_values[!is.na(directionality)]
  z_directional <- abs(sum(stats::qnorm(p_values_directional/2)*directionality[!is.na(directionality)]))
  
  # Sum the non-directional z-values
  z_nondirectional <- abs(sum(stats::qnorm(p_values[is.na(directionality)]/2)))
  
  # Combine both
  z_values <- c(z_directional, z_nondirectional)
  sum(z_values)/sqrt(k)
}


strubesMethod <- function (p_values){
  # Acquiring the unadjusted z-value from Stouffer's method
  stouffer_z <- c()
  for(i in 1:length(p_values[,1])){
    stouffer_z <- c(stouffer_z,stouffersMethod(p_values[i,]))
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


strubesDirectional <- function (p_values, scores_direction, constraints_vector){
  # Acquiring the unadjusted z-value from Stouffer's method
  stouffer_z <- c()
  for(i in 1:length(p_values[,1])){
    stouffer_z <- c(stouffer_z,stouffersDirectional(p_values[i,], scores_direction[i,],constraints_vector))
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
