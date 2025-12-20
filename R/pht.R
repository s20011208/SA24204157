#' @title Random projection independence test
#' 
#' @param x a numeric vector or matrix
#' @param y a numeric vector or matrix
#' @param num.permutations the number of permutation replications. 
#' When \code{num.permutations = 0}, the function just returns the statistic. 
#' Default: \code{num.permutations = 199}.
#'
#' @return If \code{num.permutations > 0}, the function returns a \code{htest} class object containing the following components:
#' If \code{num.permutations = 0}, the function returns a statistic value.
#' 
#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib pit, .registration = TRUE
#' 
#' @examples
#' set.seed(1)
#' num <- 30
#' 
#' ## univariate case:
#' x <- rnorm(num)
#' y <- rnorm(num)
#' pht.test(x, y)
#' y <- x + rnorm(num)
#' pht.test(x, y)
#' 
#' ## multivariate case:
#' p <- 3
#' x <- matrix(rnorm(num * p), nrow = num)
#' y <- matrix(rnorm(num * p), nrow = num)
#' pht.test(x, y)
#' y <- x + matrix(rnorm(num * p), nrow = num)
#' pht.test(x, y)
pht.test <- function(x, y, num.permutations = 99, seed = 1) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  R <- as.integer(num.permutations)
  set.seed(seed)
  res <- rcpp_pht(x, y, R)
  set.seed(NULL)
  res
}

.onUnload <- function (libpath) {
  library.dynam.unload("pit", libpath)
}


