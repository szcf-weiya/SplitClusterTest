# Wilcoxon rank sum test (adapted from `limma::rankSumTestWithCorrelation`)

Wilcoxon rank sum test (adapted from
`limma::rankSumTestWithCorrelation`)

## Usage

``` r
rankSumTestWithCorrelation(index, statistics, correlation = 0, df = Inf)
```

## Arguments

- index:

  any index vector such that `statistics[index]` contains the values of
  the statistic for the test group.

- statistics:

  numeric vector giving values of the test statistic.

- correlation:

  numeric scalar, average correlation between cases in the test group.
  Cases in the second group are assumed independent of each other and
  other the first group.

- df:

  degrees of freedom which the correlation has been estimated.

## Value

Numeric vector of length 4 containing the `left.tail` and `right.tail`
p-values and the corresponding test statistics.
