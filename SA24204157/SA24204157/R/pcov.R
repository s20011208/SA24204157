#' @title Projection correlation independence test
#' 
#' @description Permutation based projection correlation independence test
#'
#' @inheritParams pht.test
#' 
#' @param test.statistics a logical value. If \code{test.statistics = TRUE}, an adjusted test statistic would be computed. 
#' Emprically, the adjusted projection correlation can improve the statistical power but increase computational burden. 
#' Default: \code{test.statistics = FALSE}.
#' 
#' @references Liping Zhu, Kai Xu, Runze Li, Wei Zhong, Projection correlation between two random vectors, Biometrika, Volume 104, Issue 4, December 2017, Pages 829-843, https://doi.org/10.1093/biomet/asx043
#' 
#' @inherit pht.test return
#' 
#' @export
#' 
#' @examples
#' set.seed(1)
#' num <- 30
#'
#' ## univariate case:
#' x <- rnorm(num)
#' y <- rnorm(num)
#' pcov.test(x, y)
#' y <- x + rnorm(num)
#' pcov.test(x, y)
#'
#' ## multivariate case:
#' p <- 3
#' x <- matrix(rnorm(num * p), nrow = num)
#' y <- matrix(rnorm(num * p), nrow = num)
#' pcov.test(x, y)
#' y <- x + matrix(rnorm(num * p), nrow = num)
#' pcov.test(x, y)
pcov.test <- function(x, y, num.permutations = 99, test.statistics = FALSE, seed = 1) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  R <- as.integer(num.permutations)
  term_S2 <- as.integer(test.statistics)
  set.seed(seed)
  res <- rcpp_pcov(x, y, R, term_S2)
  set.seed(NULL)
  res
}

