# Wrapper for the naive double-dipping method

Wrapper for the naive double-dipping method

## Usage

``` r
dd(sce, params, test.use = "t", q = 0.05)
```

## Arguments

- sce:

  A SingleCellExperiment object

- params:

  parameters passed to the clustering step `perform_clustering`

- test.use:

  The hypothesis testing for DE test. Possible choices: `t`, `wilcox`,
  `poisson`, and `negbinom`

- q:

  nominal FDR level

## Value

a vector of selected significant features
