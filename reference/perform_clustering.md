# Perform Clustering on a SingleCellExperiment object

Perform Clustering on a SingleCellExperiment object

## Usage

``` r
perform_clustering(
  sce,
  ret_obj = FALSE,
  normalized_method = "LogNormalize",
  use.kmeans = FALSE,
  kmeans.nstart = 1,
  resolution = 0.5,
  step = 0.05,
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
  signal_measurement = "tstat",
  verbose = FALSE
)
```

## Arguments

- sce:

  A SingleCellExperiment object

- ret_obj:

  (Default: `FALSE`) whether to return the `sce` after operation instead
  of returning the signal measurement

- normalized_method:

  Normalization method. Possible choices: `LogNormalize`, `sct`, or
  `none`

- use.kmeans:

  if `TRUE`, then clustering using kmeans with two clusters, otherwise,
  use
  [`Seurat::FindClusters`](https://satijalab.org/seurat/reference/FindClusters.html)
  and need to find `resolution` to achieve two clusters

- kmeans.nstart:

  the `nstart` parameter for `kmeans`

- resolution:

  the `resolution` parameter for
  [`Seurat::FindClusters`](https://satijalab.org/seurat/reference/FindClusters.html)

- step:

  To find a resolution with desired number of clusters, the step size
  for changing the `resolution` parameter

- maxiter:

  maximum iteration for searching `resolution` to have two clusters

- louvain.nstart:

  the `n.start` parameter for
  [`Seurat::FindClusters`](https://satijalab.org/seurat/reference/FindClusters.html)

- louvain.alg:

  the `algorithm` parameter for
  [`Seurat::FindClusters`](https://satijalab.org/seurat/reference/FindClusters.html)

- seed.cluster:

  the `random.seed` parameter for
  [`Seurat::FindClusters`](https://satijalab.org/seurat/reference/FindClusters.html)

- kmeans.whiten:

  whether to whitening for the kmeans clustering

- pca.whiten:

  whether to whitening for PCA if using Louvain algorithms

- npcs:

  number of PCs, the parameter `npcs` used in
  [`Seurat::RunPCA`](https://satijalab.org/seurat/reference/RunPCA.html)

- k.param:

  the parameter `k.param` used in
  [`Seurat::FindNeighbors`](https://satijalab.org/seurat/reference/FindNeighbors.html)

- test.use:

  The hypothesis testing for DE test. Possible choices: `t`, `wilcox`,
  `poisson`, and `negbinom`

- test.from.raw:

  whether to perform the testing on the original data without
  normalization or on the normalized data

- signal_measurement:

  the signal measurement. Possible choices: the test statistic `tstat`,
  or the signed p-value `pval`

- verbose:

  whether to print internal log messages

## Value

a vector of signal measurements for each feature
