# Conduct multiple data splitting in parallel

Step 1 (in parallel) of `mds`

## Usage

``` r
mds1_parallel(sce, M = 2, ncores = 10, ...)
```

## Arguments

- sce:

  A SingleCellExperiment object

- M:

  The number of data splitting

- ncores:

  The number of cores used for parallel computing. It is better to be a
  factor of `M`.

- ...:

  Arguments passed to `ds`

## Value

a list of the mirror statistics from multiple data splitting
