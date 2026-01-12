# MDS procedure for matrix

MDS procedure for matrix

## Usage

``` r
# S3 method for class 'matrix'
mds(
  x,
  M = 10,
  q = 0.05,
  ncores = 1,
  tied.method = "fair",
  sum_in_denom = TRUE,
  ...
)
```

## Arguments

- x:

  data matrix

- M:

  number of repetitions

- q:

  nominal FDR level

- ncores:

  number of cores

- tied.method:

  handling the tied values when determing the cutoff

- sum_in_denom:

  way for calculating the inclusion rate

- ...:

  optional arguments passed to `ds.matrix`

## Value

a list of mirror statistics
