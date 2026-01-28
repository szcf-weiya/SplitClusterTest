#' DS procedure
#' @param x data matrix or SingleCellExperiment object
#' @param ... optional arguments
#' @export
ds = function(x, ...) UseMethod("ds")

#' MDS procedure
#' @param x data matrix or SingleCellExperiment object
#' @param ... optional arguments
#' @export
mds = function(x, ...) UseMethod("mds")

#' DS procedure for Matrix
#' @param x data matrix
#' @param kmeans_whiten whether to perform the whitening procedure for kmean
#' @param Sigma the covariance matrix, if null, estimate from the data
#' @param q the nominal FDR level
#' @param debias whether to debias (useful when correlation is high)
#' @return a list with the selected set, the mirror statistics, two signal measurements for two halves
#' @export
ds.matrix = function(x, kmeans_whiten = FALSE, Sigma = NULL, q = 0.05, debias = FALSE) {
  n = nrow(x)
  n2 = floor(n/2)
  idx1 = sample(1:n, n2, replace = FALSE)
  idx2 = setdiff(1:n, idx1)
  x1 = x[idx1, ]
  x2 = x[idx2, ]
  if (kmeans_whiten) {
    p = ncol(x)
    if (is.null(Sigma)) {
      Sigma = est.Sigma(x) + 1e-6 * diag(p)
    }
    ev = eigen(Sigma)
    idx = ev$values > 0
    xc = x %*% ev$vectors[, idx] %*% diag(1 / sqrt(ev$values[idx])) %*% t(ev$vectors[, idx])
    x1c = xc[idx1, ]
    x2c = xc[idx2, ]
  } else {
    x1c = x1
    x2c = x2
  }
  cl1 = kmeans(x1c, 2)$cluster
  cl2 = kmeans(x2c, 2)$cluster
  d1 = cluster_diff(x1, cl1)
  d2 = cluster_diff(x2, cl2)
  if (debias) {
    d1 = debias_symmetry(d1)
    d2 = debias_symmetry(d2)
  }
  ms = mirror_stat(d1, d2)
  tau = calc_tau(ms, q)
  m_sel = which(ms >= tau) # >= consistent with the FDR calculation
  names(ms) = 1:length(ms) # to be compatible with gene names
  list(sel_set = m_sel, ms = ms, d1 = d1, d2 = d2)
}


#' DS procedure for Seurat object
#' @param sce A SingleCellExperiment object
#' @param params1 list of parameters to be passed to `perform_clustering` for the first half of the data splitting
#' @param params2 list of parameters to be passed to `perform_clustering` for the second half of the data splitting.
#' @param seed random seed for reproducibility
#' @inheritParams perform_clustering
#' @return a vector of mirror statistics
#' @importFrom stats runif
#' @import Seurat
#' @export
ds.SingleCellExperiment = function(sce, params1 = NULL, params2 = NULL,
              seed = NA,
              test.use = "t",
              signal_measurement = "tstat",
              test.from.raw = TRUE) {
  if (is.na(seed)) {
    seed = round(runif(1) * 1000)
    cat("seed = ", seed, "\n")
  }
  set.seed(seed)
  n = ncol(sce)
  idx1 = sample(1:n, round(n/2), replace = FALSE)
  idx2 = setdiff(1:n, idx1)
  sce1 = sce[, idx1]
  sce2 = sce[, idx2]
  ## common arguments for each half
  params1$sce = sce1
  params2$sce = sce2
  params1$test.use = test.use
  params2$test.use = test.use
  params1$signal_measurement = signal_measurement
  params2$signal_measurement = signal_measurement
  params1$test.from.raw = test.from.raw
  params2$test.from.raw = test.from.raw
  # just take different random seed
  params1$seed.cluster = seed + 10000
  params2$seed.cluster = seed + 20000
  cat("\n----- data splitting (1st half) -----\n")
  signed_logpval1.raw = do.call(perform_clustering, params1)
  cat("\n----- data splitting (2nd half) ----\n")
  signed_logpval2.raw = do.call(perform_clustering, params2)

  all_names = union(names(signed_logpval1.raw), names(signed_logpval2.raw))
  signed_logpval1 = structure(numeric(length(all_names)), names = all_names)
  signed_logpval2 = signed_logpval1
  signed_logpval1[names(signed_logpval1.raw)] = signed_logpval1.raw
  signed_logpval2[names(signed_logpval2.raw)] = signed_logpval2.raw
  signed_logpval1[is.na(signed_logpval1)] = 0
  signed_logpval2[is.na(signed_logpval2)] = 0
  lbl.sgn = sum(signed_logpval1 * signed_logpval2, na.rm = T)
  m = (abs(signed_logpval1) + abs(signed_logpval2)) * sign(signed_logpval1) * sign(signed_logpval2) * sign(lbl.sgn)
  m
}

#' Wrapper for the naive double-dipping method
#' @param sce A SingleCellExperiment object
#' @param params parameters passed to the clustering step `perform_clustering`
#' @inheritParams perform_clustering
#' @param q nominal FDR level
#' @return a vector of selected significant features
#' @importFrom stats p.adjust
#' @export
dd = function(sce, params, test.use = "t", q = 0.05) {
  params$sce = sce
  params$ret_obj = TRUE
  sce = do.call(perform_clustering, params)
  res = tryCatch({FindMarkers(sce,
                    ident.1 = 0,
                    ident.2 = 1,
                    min.pct = 0,
                    logfc.threshold = 0, test.use = test.use)},
                 error = function(e) e)
  if (!isa(res, "try-error"))
    rownames(res)[p.adjust(res$p_val, "BH") < q]
  else {
    cat("failed to clustering into two groups\n")
  }
}

#' Calculate the cutoff the mirror statistics
#' @param ms Array of mirror statistics
#' @param q nominal FDR level (default: 0.05)
#' @return a cutoff
#' @export
calc_tau = function(ms, q = 0.05) {
  ts = unique(sort(abs(ms)))
  cutoff_set = max(ms) + 1 # since the selection is >=, avoid max(ms) to be selected if no proper cutoff
  for (t in ts) {
    curr_fdr = (1 + sum(ms <= -t)) / max(1, sum(ms >= t))
    if (curr_fdr <= q)
      cutoff_set = c(cutoff_set, t)
  }
  min(cutoff_set)
}

#' Conduct multiple data splitting
#'
#' Step 1 of `mds`
#'
#' @param sce A SingleCellExperiment object
#' @param M The number of data splitting
#' @param ... Arguments passed to `ds`
#' @return a list of the mirror statistics from multiple data splitting
#' @export
mds1 = function(sce, M = 2, ...) {
  mss = list()
  for (i in 1:M) {
    cat("\n===== Multiple Data Splitting: ", i, "/", M, " =====\n")
    mss[[i]] = ds(sce, seed = i, ...)
  }
  return(mss)
}

#' Conduct multiple data splitting in parallel
#'
#' Step 1 (in parallel) of `mds`
#'
#' @param sce A SingleCellExperiment object
#' @param M The number of data splitting
#' @param ncores The number of cores used for parallel computing. It is better to be a factor of `M`.
#' @param ... Arguments passed to `ds`
#' @return a list of the mirror statistics from multiple data splitting
#' @export
mds1_parallel <- function(sce, M = 2, ncores = 10, ...) {
  mss <- pbmcapply::pbmclapply(1:M, function(i) {
    ds(sce, seed = i, ...)
  }, mc.cores = ncores)
  return(mss)
}

#' Calculate the inclusion rate
#' @param mss a list of the selection of multiple data splitting
#' @param q the nominal FDR level
#' @param sum_in_denom use the modified inclusion rate (`TRUE`) or the original inclusion rate (`FALSE`)
#' @return the inclusion rate for each feature
calc_inc_rate = function(mss, q = 0.05, sum_in_denom = TRUE) {
  M = length(mss)
  # align names
  all_names = unique(unlist(lapply(mss, names)))
  taus = numeric(M)
  p = length(all_names)
  rs = matrix(0, nrow = M, ncol = p)
  colnames(rs) = all_names
  for (i in 1:M) {
    if (length(mss[[i]]) > 0) {
      taus[i] = calc_tau(mss[[i]], q = q)
      rs[i, names(mss[[i]])] = mss[[i]] >= taus[i]
      if (!sum_in_denom) rs[i, ] = rs[i, ] / max(1, sum(rs[i, ])) # rs[i, ] can be zero
    }
  }
  rs[is.na(rs)] = 0
  rss = rowSums(rs)
  r_denom = sum(rss)
  r_num = colSums(rs)
  if (sum_in_denom) {
    # intuitively, the selected items should appear in most of splits
    if (max(r_num) <= M * 0.5)
      return(numeric(p))
  }
  if (sum_in_denom)
    inc_rate = r_num / r_denom
  else
    inc_rate = r_num / M
  nzero = sum(rss == 0)
  if (nzero >= M * 0.5 )
    return(numeric(p))
  return(inc_rate)
}

#' Perform selection based on inclusion rate
#' @param inc_rate inclusion rate
#' @param q nominal FDR level
#' @param tied.method the method for handling the tied values. Possible choices: \describe{
#'  \item{`include`}{include the tied rate}
#'  \item{`exclude`}{exclude the tied rate}
#'  \item{`fair` (default)}{a compromise way which compare the number of tied values before and after the cutoff}
#' }
#' @param tol tolerance when comparing the equality of two doubles
#' @return the selected relevant features
sel_inc_rate = function(inc_rate, q = 0.05, tied.method = "fair", tol = 1e-10) {
  p = length(inc_rate)
  sort_inc_rate = sort(inc_rate)
  curr_rate = 0
  for (i in 1:(p+1)) {
    if (curr_rate > q) {
      cat("use the tie.method = ", tied.method, "\n")
      if (tied.method == "include") { # always include the tie value, it can be powerful, but fdr might not good
        cat("include the tied cutoff inc_rate = ", sort_inc_rate[i-1], "\n")
        return(which(inc_rate >= sort_inc_rate[i-1]))
      } else if (tied.method == "exclude") {
        cat("cutoff inc_rate = ", sort_inc_rate[i-1], "\n")
        if (abs(sort_inc_rate[i-1] - sort_inc_rate[p]) < tol) # only when the inc rate is large
          return(which(inc_rate >= sort_inc_rate[i-1]))
        else
          return(which(inc_rate > sort_inc_rate[i-1]))
      } else {
        if (i <= p)
          ntie = sum(abs(sort_inc_rate[i:p] - sort_inc_rate[i-1]) < tol)
        else
          ntie = 0 # if there are many tied values, the for loop should be terminated before i = p+1. Now if it only large than q*r_denom at the last point, return empty.
        ntie0 = sum(abs(sort_inc_rate[1:(i-1)] - sort_inc_rate[i-1]) < tol)
        if (ntie >= ntie0) {
          return(which(inc_rate >= sort_inc_rate[i-1]))
        } else {
          return(which(inc_rate > sort_inc_rate[i-1]))
        }
      }
    } else {
      if (i == p+1) { # the total summation did not achieve the rate
        return(which(inc_rate > 0))
      }
      curr_rate = curr_rate + sort_inc_rate[i]
    }
  }
}

#' MDS procedure for matrix
#' @param x data matrix
#' @param M number of repetitions
#' @param q nominal FDR level
#' @param ncores number of cores
#' @param tied.method handling the tied values when determing the cutoff
#' @param sum_in_denom way for calculating the inclusion rate
#' @param ... optional arguments passed to `ds.matrix`
#' @return a list of mirror statistics
#' @export
mds.matrix = function(x, M = 10, q = 0.05,
                      ncores = 1,
                      tied.method = "fair", sum_in_denom = TRUE, ...) {
  mss = list()
  if (ncores == 1) {
    for (i in 1:M) {
      mss[[i]] = ds.matrix(x, ...)$ms
    }
  } else {
    mss <- pbmcapply::pbmclapply(1:M, function(i) {
      ds(x, ...)$ms
    }, mc.cores = ncores)
  }
  sel = mds2(mss, q = q, tied.method = tied.method, sum_in_denom = sum_in_denom)
  return(sel)
}

#' Aggregate multiple data splitting results, and return a selection set
#' @param mss A list, which contains the selections from multiple data splitting procedures
#' @param q The nominal FDR level (default: 0.05)
#' @param tied.method The method for handling the tied values.
#' @param sum_in_denom when calculating the inclusion rate, use the modified inclusion rate (`TRUE`) or the original inclusion rate (`FALSE`)
#' @return the selected relevant features
#' @export
mds2 = function(mss, q = 0.05, tied.method = "fair", sum_in_denom = TRUE) {
  mss = mss[sapply(mss, function(x) !isa(x, "try-error"))] # possible due to mds_parallel
  M = length(mss)
  if (M < 1) {
    cat("!!!!!!WARNING: no valid ds split!!!!\n\n")
    return(character(0))
  }
  inc_rate = calc_inc_rate(mss, q = q, sum_in_denom = sum_in_denom)
  ret = sel_inc_rate(inc_rate, q = q, tied.method = tied.method)
  return(names(ret))
}

#' Multiple data splitting
#' @param sce A SingleCellExperiment object
#' @param ncores If larger than 1, then compute in parallel with `ncores` cores.
#' @param M The number of data splitting
#' @inheritParams mds2
#' @param ... arguments passed to `mds1` (or `mds1_parallel`)
#' @return the selected relevant features
#' @export
mds.SingleCellExperiment = function(sce, ncores = 1, M = 10, q = 0.05,
               tied.method = "fair", sum_in_denom = TRUE, ...) {
  if (ncores > 1) {
    mss = mds1_parallel(sce, M = M, ncores = ncores, ...)
  } else {
    mss = mds1(sce, M = M, ...)
  }
  sel = mds2(mss, q = q, tied.method = tied.method, sum_in_denom = sum_in_denom)
  return(sel)
}

#' Perform Clustering on a SingleCellExperiment object
#' @param sce A SingleCellExperiment object
#' @param ret_obj (Default: `FALSE`) whether to return the `sce` after operation instead of returning the signal measurement
#' @param normalized_method Normalization method. Possible choices: `LogNormalize`, `sct`, or `none`
#' @param use.kmeans if `TRUE`, then clustering using kmeans with two clusters, otherwise, use `Seurat::FindClusters` and need to find `resolution` to achieve two clusters
#' @param kmeans.nstart the `nstart` parameter for `kmeans`
#' @param resolution the `resolution` parameter for `Seurat::FindClusters`
#' @param step To find a resolution with desired number of clusters, the step size for changing the `resolution` parameter
#' @param maxiter maximum iteration for searching `resolution` to have two clusters
#' @param louvain.nstart the `n.start` parameter for `Seurat::FindClusters`
#' @param louvain.alg the `algorithm` parameter for `Seurat::FindClusters`
#' @param seed.cluster the `random.seed` parameter for `Seurat::FindClusters`
#' @param kmeans.whiten whether to whitening for the kmeans clustering
#' @param pca.whiten whether to whitening for PCA if using Louvain algorithms
#' @param npcs number of PCs, the parameter `npcs` used in `Seurat::RunPCA`
#' @param k.param the parameter `k.param` used in `Seurat::FindNeighbors`
#' @param test.use The hypothesis testing for DE test. Possible choices: `t`, `wilcox`, `poisson`, and `negbinom`
#' @param test.from.raw whether to perform the testing on the original data without normalization or on the normalized data
#' @param signal_measurement the signal measurement. Possible choices: the test statistic `tstat`, or the signed p-value `pval`
#' @param verbose whether to print internal log messages
#' @return a vector of signal measurements for each feature
#' @importFrom stats cov kmeans
#' @export
perform_clustering = function(sce,
                              ret_obj = FALSE,
                              normalized_method = "LogNormalize",
                              use.kmeans = FALSE, kmeans.nstart = 1,
                              resolution = 0.5, step = 0.05,
                              maxiter = 40,
                              louvain.nstart = 1,
                              louvain.alg = 1,
                              seed.cluster = 0,
                              kmeans.whiten = FALSE,
                              pca.whiten = FALSE,
                              npcs = 10,
                              k.param = 20,
                              test.use = "t",
                              test.from.raw = FALSE,
                              signal_measurement = "tstat", # or signed pvalue
                              verbose = FALSE
) {
  obj = CreateSeuratObject(counts = SingleCellExperiment::counts(sce))
  if (normalized_method == "sct") {
    obj <- suppressWarnings({SCTransform(obj, verbose = verbose, return.only.var.genes = FALSE)}) ## suppress "Warning: The following arguments are not used: norm.method"
  } else if (normalized_method == "LogNormalize") {
    obj <- NormalizeData(object = obj, normalization.method = normalized_method, verbose = verbose)
    obj <- FindVariableFeatures(object = obj, verbose = verbose)
    obj <- ScaleData(object = obj, verbose = verbose)
  } else {
    # do nothing for the normalization
    obj <- FindVariableFeatures(object = obj, verbose = verbose, nfeatures = nrow(obj))
    obj <- ScaleData(object = obj, verbose = verbose)
  }
  # Centering and scaling data matrix
  n = ncol(obj)
  cl = sample(1:2, n, replace = T) # dummy
  ncl = 2
  tryCatch({
    obj <- RunPCA(object = obj, verbose = verbose, weight.by.var = !pca.whiten, npcs = npcs)
    obj <- FindNeighbors(object = obj, verbose = verbose, dims = 1:npcs, k.param = k.param)
    obj <- FindClusters(object = obj, resolution = resolution, verbose = verbose,
                        random.seed = seed.cluster, n.start = louvain.nstart,
                        algorithm = louvain.alg)
    cl = obj@meta.data$seurat_clusters
    ncl = length(unique(cl))
  }, error = function(e) {
    print(e)
    cat("use k-means instead\n")
    use.kmeans <<- TRUE
  })
  visited.resolutions = c(resolution)
  visited.ncl = c(ncl)
  iter = 0
  resolution.left = NA
  resolution.right = NA
  while (TRUE) {
    nvisit = length(visited.ncl)
    iter = iter + 1
    if (iter > maxiter || use.kmeans) {
      if (use.kmeans) {
        cat("Use kmeans to force two clusters\n")
      } else {
        cat("cannot cluster into 2 clusters using FindClusters\n")
      }
      set.seed(seed.cluster)
      if (kmeans.whiten) {
        x = t(obj[[obj@active.assay]]@scale.data)
        Sigma = cov(x)
        ev = eigen(Sigma)
        eval = ev$values
        eval[eval <= 0] = 1e-16
        x2 = x %*% ev$vectors %*% diag(1 / sqrt(eval))
        km_res = kmeans(x2, centers = 2, nstart = kmeans.nstart)
      }
      else {
        km_res <- kmeans(t(obj[[obj@active.assay]]@scale.data), centers = 2, nstart = kmeans.nstart)
      }
      # NOTE: pam returns `clustering` while kmeans returns `cluster`. We can simply use `cluster` as the vector name.
      if (min(table(km_res$cluster)) < 3) # avoid error "Cell group 2 has fewer than 3 cells"
      {
        cat("\nrandomly set clustering label due to 'One group has fewer than 3 cells'\n")
        n = length(km_res$cluster)
        dummy_cl = sample(1:2, n, replace = T)
        km_res$cluster[1:n] = dummy_cl # keep the names
      }
      obj@meta.data$seurat_clusters <- factor(km_res$cluster - 1)
      names(obj@meta.data$seurat_clusters) <- colnames(obj)
      obj@active.ident <- factor(km_res$cluster - 1)
      cl = obj@meta.data$seurat_clusters
      ncl = length(unique(cl))
      break
    }
    if (ncl == 2) break
    else if (ncl < 2) { # increase resolution
      resolution.left = resolution
      resolution = resolution + step
    } else {
      resolution.right = resolution
      resolution = max(resolution - step, resolution / 2)
    }
    if (is.na(resolution.left) & is.na(resolution.right)) {
      resolution = (resolution.left + resolution.right) / 2
      next
    }
    obj <- FindClusters(object = obj, resolution = resolution, verbose = verbose,
                        random.seed = seed.cluster, n.start = louvain.nstart,
                        algorithm = louvain.alg)
    cl = obj@meta.data$seurat_clusters
    ncl = length(unique(cl))
    visited.resolutions = c(visited.resolutions, resolution)
    visited.ncl = c(visited.ncl, ncl)
    if (verbose) {
      cat("try resolution = ", resolution,  " step = ", step, " ")
      cat("ncl = ", ncl, "\n")
    }
  }
  if (ret_obj) return(obj)
  if (ncl == 1) {
    pval = numeric(nrow(obj))
    names(pval) = rownames(obj)
    return(pval)
  }
  if (test.from.raw) {
    cell.names = colnames(obj[["RNA"]])
    markers = myFindMarkers(obj[["RNA"]],
                            cells.1 = cell.names[obj@meta.data$seurat_clusters == 0],
                            cells.2 = cell.names[obj@meta.data$seurat_clusters == 1],
                            min.pct = 0,
                            logfc.threshold = 0,
                            test.use = test.use)
  } else {
    markers = myFindMarkers(obj,
                            ident.1 = 0,
                            ident.2 = 1,
                            min.pct = 0,
                            logfc.threshold = 0,
                            test.use = test.use)
  }
  if (signal_measurement == "tstat")
    pval = markers$tstat
  else
    pval = -log10(markers$p_val) * sign(markers$tstat)
  names(pval) = rownames(markers)
  pval
}
