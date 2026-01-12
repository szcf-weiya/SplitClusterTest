# DS procedure for Seurat object

DS procedure for Seurat object

## Usage

``` r
# S3 method for class 'SingleCellExperiment'
ds(
  sce,
  params1 = NULL,
  params2 = NULL,
  seed = NA,
  test.use = "t",
  signal_measurement = "tstat",
  test.from.raw = TRUE
)
```

## Arguments

- sce:

  A SingleCellExperiment object

- params1:

  list of parameters to be passed to `perform_clustering` for the first
  half of the data splitting

- params2:

  list of parameters to be passed to `perform_clustering` for the second
  half of the data splitting.

- seed:

  random seed for reproducibility

- test.use:

  The hypothesis testing for DE test. Possible choices: `t`, `wilcox`,
  `poisson`, and `negbinom`

- signal_measurement:

  the signal measurement. Possible choices: the test statistic `tstat`,
  or the signed p-value `pval`

- test.from.raw:

  whether to perform the testing on the original data without
  normalization or on the normalized data

## Value

a vector of mirror statistics
