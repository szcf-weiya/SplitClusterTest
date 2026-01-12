# Aggregate multiple data splitting results, and return a selection set

Aggregate multiple data splitting results, and return a selection set

## Usage

``` r
mds2(mss, q = 0.05, tied.method = "fair", sum_in_denom = TRUE)
```

## Arguments

- mss:

  A list, which contains the selections from multiple data splitting
  procedures

- q:

  The nominal FDR level (default: 0.05)

- tied.method:

  The method for handling the tied values.

- sum_in_denom:

  when calculating the inclusion rate, use the modified inclusion rate
  (`TRUE`) or the original inclusion rate (`FALSE`)

## Value

the selected relevant features
