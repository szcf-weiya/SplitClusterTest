# Generate simulation data with two gaussians

Generate simulation data with two gaussians

## Usage

``` r
gen_data_normal(
  n = 100,
  p = 200,
  delta = 3,
  prop = 0.1,
  rho = 0,
  prop_cl1 = 0.5
)
```

## Arguments

- n:

  sample size

- p:

  number of features

- delta:

  signal strength

- prop:

  proportion of non-null features

- rho:

  correlation parameter in the covariance matrix

- prop_cl1:

  proportion of the first cluster
