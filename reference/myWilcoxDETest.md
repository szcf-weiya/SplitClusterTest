# Modified Seurat::WilcoxDETest by extending the output with statistics in addition to p-values

Modified Seurat::WilcoxDETest by extending the output with statistics in
addition to p-values

## Usage

``` r
myWilcoxDETest(data.use, cells.1, cells.2, verbose = TRUE, ...)
```

## Arguments

- data.use:

  Data matrix to test

- cells.1:

  Group 1 cells

- cells.2:

  Group 2 cells

- verbose:

  Print a progress bar

- ...:

  Extra parameters passed to wilcox.test

## Value

Returns a p-value ranked matrix of putative differentially expressed
features
