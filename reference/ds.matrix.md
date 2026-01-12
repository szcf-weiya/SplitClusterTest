# DS procedure for Matrix

DS procedure for Matrix

## Usage

``` r
# S3 method for class 'matrix'
ds(x, kmeans_whiten = FALSE, Sigma = NULL, q = 0.05, debias = FALSE)
```

## Arguments

- x:

  data matrix

- kmeans_whiten:

  whether to perform the whitening procedure for kmean

- Sigma:

  the covariance matrix, if null, estimate from the data

- q:

  the nominal FDR level

- debias:

  whether to debias (useful when correlation is high)

## Value

a list with the selected set, the mirror statistics, two signal
measurements for two halves
