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
  L = c(rep(1, n1), rep(2, n - n1))
  return(list(X = x, L = L))
}

#' Generate Poisson data with latent structure
#'
#' @param n Number of observations
#' @param p Number of features
#' @param delta Effect size parameter
#' @param prop_imp Proportion of important features (default: 0.1)
#' @param rho AR(1) correlation parameter (default: 0.5)
#' @param block_size Size of correlation blocks (default: 10)
#' @param type Type of latent variable: "discrete" (binary) or "continuous" (default: "discrete")
#' @param sigma Additional noise standard deviation (default: 0)
#' @importFrom stats rnorm
#' @return List containing X (data matrix) and L (latent variable)
#' @export
gen_data_pois <- function(n, p, delta = 3, prop_imp = 0.1, rho = 0.5,
                          block_size = 10, type = "discrete", sigma = 0) {

  # Generate latent variable L
  if (type == "discrete") {
    L <- sample(0:1, n, replace = TRUE)
  } else {
    L <- rnorm(n)
    L <- L - mean(L)  # Center the latent variable
  }

  # Create factor loadings Fs
  p1 <- as.integer(p * prop_imp)
  p0 <- p - p1
  Fs <- c(rep(delta, p1), rep(0, p0))

  # Create L * Fs' (outer product)
  LFT <- outer(L, Fs)

  # Create logLambda with intercept log(3) ~ 1.098612
  logLambda <- LFT + log(3) + matrix(rnorm(n * p, sd = sigma), nrow = n, ncol = p)

  # Convert to Lambda
  Lambda <- exp(logLambda)

  # Generate Poisson data using the matrix version
  X <- gen_data_pois.matrix(Lambda, rho = rho, block_size = block_size)

  return(list(X = X, L = L))
}

#' Matrix method for gen_data_pois (S3)
#' @param Lambda parameter matrix
#' @param rho correlation parameter
#' @param block_size the copula correlation is set up in each block (mainly for speed up copula for multivariate distribution)
#' @return data matrix of the same size of `Lambda`
#' @importFrom copula normalCopula P2p rMvdc mvdc
gen_data_pois.matrix <- function(Lambda, rho = 0.5, block_size = 10) {
  n <- nrow(Lambda)
  p <- ncol(Lambda)

  # Check if p is divisible by block_size
  if (p %% block_size != 0) {
    stop("p must be divisible by block_size")
  }

  Sigma <- create_ar1_matrix(rho, block_size)
  pp <- as.integer(p / block_size)
  X <- matrix(0, nrow = n, ncol = p)

  cop <- normalCopula(P2p(Sigma), dim = block_size, dispstr = "un")
  # Use Gaussian copula approach
  for (i in 1:n) {
    for (j in 1:pp) {
      paras = vector("list", length = block_size)
      for (k in 1:block_size) {
        paras[[k]] = list(lambda = Lambda[i, (j-1)*block_size + k])
      }
      X[i, ((j-1)*block_size+1) : (j*block_size)] = rMvdc(1, mvdc(cop, rep("pois", block_size), paras))
    }
  }
  return(X)
}

# Create AR(1) correlation matrix
create_ar1_matrix <- function(rho, size) {
  Sigma <- matrix(0, nrow = size, ncol = size)
  for (i in 1:size) {
    for (j in 1:size) {
      Sigma[i, j] <- rho^(abs(i - j))
    }
  }
  return(Sigma)
}
