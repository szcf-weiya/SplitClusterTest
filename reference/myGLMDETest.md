# Modified Seurat::GLMDETest by extending the output with statistics in addition to p-values

Modified Seurat::GLMDETest by extending the output with statistics in
addition to p-values

## Usage

``` r
myGLMDETest(
  data.use,
  cells.1,
  cells.2,
  min.cells = 3,
  latent.vars = NULL,
  test.use = NULL,
  verbose = TRUE
)
```

## Arguments

- data.use:

  Data to test

- cells.1:

  Group 1 cells

- cells.2:

  Group 2 cells

- min.cells:

  Minimum number of cells threshold

- latent.vars:

  Latent variables to test

- test.use:

  parameterizes the glm

- verbose:

  Print progress bar

## Value

Returns a p-value ranked matrix of putative differentially expressed
