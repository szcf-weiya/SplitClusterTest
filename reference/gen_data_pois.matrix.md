# Matrix method for gen_data_pois (S3)

Matrix method for gen_data_pois (S3)

## Usage

``` r
gen_data_pois.matrix(Lambda, rho = 0.5, block_size = 10)
```

## Arguments

- Lambda:

  parameter matrix

- rho:

  correlation parameter

- block_size:

  the copula correlation is set up in each block (mainly for speed up
  copula for multivariate distribution)

## Value

data matrix of the same size of `Lambda`
