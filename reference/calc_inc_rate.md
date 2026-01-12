# Calculate the inclusion rate

Calculate the inclusion rate

## Usage

``` r
calc_inc_rate(mss, q = 0.05, sum_in_denom = TRUE)
```

## Arguments

- mss:

  a list of the selection of multiple data splitting

- q:

  the nominal FDR level

- sum_in_denom:

  use the modified inclusion rate (`TRUE`) or the original inclusion
  rate (`FALSE`)

## Value

the inclusion rate for each feature
