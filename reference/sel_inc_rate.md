# Perform selection based on inclusion rate

Perform selection based on inclusion rate

## Usage

``` r
sel_inc_rate(inc_rate, q = 0.05, tied.method = "fair", tol = 1e-10)
```

## Arguments

- inc_rate:

  inclusion rate

- q:

  nominal FDR level

- tied.method:

  the method for handling the tied values. Possible choices:

  `include`

  :   include the tied rate

  `exclude`

  :   exclude the tied rate

  `fair` (default)

  :   a compromise way which compare the number of tied values before
      and after the cutoff

- tol:

  tolerance when comparing the equality of two doubles

## Value

the selected relevant features
