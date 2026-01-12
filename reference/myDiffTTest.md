# Modified Seurat::DiffTTest by extending the output with statistics in addition to p-values

Differential expression testing using Student's t-test

## Usage

``` r
myDiffTTest(data.use, cells.1, cells.2, verbose = TRUE)
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

## Details

Identify differentially expressed genes between two groups of cells
using the Student's t-test
