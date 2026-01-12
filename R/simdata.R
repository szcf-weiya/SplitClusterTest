#' Generate simulation data with two gaussians
#' @param n sample size
#' @param p number of features
#' @param delta signal strength
#' @param prop proportion of non-null features
#' @param rho correlation parameter in the covariance matrix
#' @param prop_cl1 proportion of the first cluster
#' @importFrom MASS mvrnorm
#' @export
gen_data_normal = function(n = 100, p = 200, delta = 3, prop = 0.1, rho = 0, prop_cl1 = 0.5) {
  p1 = floor(p * prop)
  Sigma = matrix(rho, nrow = p, ncol = p)
  Sigma = Sigma + (1 - rho) * diag(p)
  x = mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  n1 = floor(n * prop_cl1)
  x[1:n1, 1:p1] = x[1:n1, 1:p1] + delta
  cl = c(rep(1, n1), rep(2, n - n1))
  return(list(x = x, cl = cl))
}
