# Generate Poisson data with latent structure

Generate Poisson data with latent structure

## Usage

``` r
gen_data_pois(
  n,
  p,
  delta = 3,
  prop_imp = 0.1,
  rho = 0.5,
  block_size = 10,
  type = "discrete",
  sigma = 0
)
```

## Arguments

- n:

  Number of observations

- p:

  Number of features

- delta:

  Effect size parameter

- prop_imp:

  Proportion of important features (default: 0.1)

- rho:

  AR(1) correlation parameter (default: 0.5)

- block_size:

  Size of correlation blocks (default: 10)

- type:

  Type of latent variable: "discrete" (binary) or "continuous" (default:
  "discrete")

- sigma:

  Additional noise standard deviation (default: 0)

## Value

List containing X (data matrix) and L (latent variable)
