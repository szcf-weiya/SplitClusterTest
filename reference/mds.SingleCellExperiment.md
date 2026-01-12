# Multiple data splitting

Multiple data splitting

## Usage

``` r
# S3 method for class 'SingleCellExperiment'
mds(
  sce,
  ncores = 1,
  M = 10,
  q = 0.05,
  tied.method = "fair",
  sum_in_denom = TRUE,
  ...
)
```

## Arguments

- sce:

  A SingleCellExperiment object

- ncores:

  If larger than 1, then compute in parallel with `ncores` cores.

- M:

  The number of data splitting

- q:

  The nominal FDR level (default: 0.05)

- tied.method:

  The method for handling the tied values.

- sum_in_denom:

  when calculating the inclusion rate, use the modified inclusion rate
  (`TRUE`) or the original inclusion rate (`FALSE`)

- ...:

  arguments passed to `mds1` (or `mds1_parallel`)

## Value

the selected relevant features
