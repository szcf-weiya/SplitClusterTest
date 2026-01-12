#' estimate the covariance matrix (assuming there are two clusters)
#' @param x data matrix
#' @return estimate covariance matrix
est.Sigma = function(x) {
  n = nrow(x)
  cl = kmeans(x, 2, nstart = 10)$cluster
  Sigma1 = cov(x[cl == 1, ])
  Sigma2 = cov(x[cl == 2, ])
  Sigma = (Sigma1 * sum(cl == 1) + Sigma2 * sum(cl == 2)) / (n - 2)
  Sigma
}

#' Calculate the difference across two clusters for each feature
#' @param x the data matrix
#' @param cl the clustering label
#' @return a vector of differences
cluster_diff = function(x, cl) {
  x1 = x[cl == 1, ]
  x2 = x[cl == 2, ]
  p = ncol(x)
  diffs = numeric(p)
  for (i in 1:p) {
    diffs[i] = t.test(x1[, i], x2[, i])$statistic
  }
  diffs
}

#' Calculate the mirror statistics
#' @param d1 signal measurement from the first half
#' @param d2 signal measurement from the second half
#' @return a vector of mirror statistics
mirror_stat = function(d1, d2) {
  (abs(d1) + abs(d2)) * sign(d1) * sign(d2) * sign(sum(d1 * d2))
}

# particularly for correlated case
# take the t-stat as an example
# assume at least prop are null features
#' Debias the statistics under the null for symmetry
#' @param t the statistics to be debiased
#' @param prop the proportion used to estimate the mean/median of statistics among the null features
#' @return the debiased statistics
#' @importFrom stats median
debias_symmetry = function(t, prop = 0.5) {
  p = length(t)
  p0 = floor(p* prop)
  idx = order(abs(t))
  b = median(t[idx][1:p0])
  t - b
}

#' Calculate the accuracy (FDR, Power, F1 score) if truth is available
#' @param est the selected relavent features
#' @param truth A vector of the true relevant features
#' @return a vector of size three: fdr, power, and F1 score
#' @export
calc_acc = function(est, truth) {
  if (length(truth) == 0) { # truth is empty
    if (length(est) == 0) {
      fdr = 0
      power = 1
      f1 = 1
    } else {
      fdr = 1
      power = 0
      f1 = 0
    }
  } else if (length(est) == 0) { # truth is not empty but return nothing
    fdr = 0
    power = 0
    f1 = 0
  } else {
    tp = length(intersect(est, truth))
    precision = tp / length(est)
    recall = tp / length(truth)
    f1 = 2 * precision * recall / (precision + recall)
    if (is.nan(f1)) f1 = 0 # 0/0
    power = recall
    fdr = 1 - precision
  }
  return(c(fdr = fdr, power = power, f1 = f1))
}

