# Two Gaussians

``` r
library(SplitClusterTest)
```

Imagine you have data from two distinct clusters, each with 200 features
and 50 observations. The interesting part: only 10% of these features
(20 features) actually differ between the groups. The remaining features
are identically distributed. The magnitude of these differences is
controlled by a parameter called `delta`.

## Generating Example Data

Let’s create such a dataset using the `gen_data_normal` function:

``` r
data = gen_data_normal(n = 100, p = 200, prop = 0.1, delta = 3)
x = data$x
cl = data$cl
```

## Identifying Non-null Features by DS

Now, let’s use our DS procedure to automatically identify which features
distinguish the two clusters:

``` r
set.seed(1)
res = ds(x)
names(res$sel_set)
#> NULL
```

Measure the accuracy of the selected set, where the true Non-null
features are `1:20`:

``` r
calc_acc(res$sel_set, 1:20)
#>   fdr power    f1 
#>     0     1     1
```

The histogram of the mirror statistics is

``` r
hist(res$ms)
```

![](two-gaussians_files/figure-html/unnamed-chunk-5-1.png)

## Identifying Non-null Features by MDS

Alternatively, to be more robust, we can apply the MDS procedure

``` r
set.seed(1)
res = mds(x, M = 10)
#> use the tie.method =  fair
```

Check the accuracy of the selected set

``` r
calc_acc(res, 1:20)
#>   fdr power    f1 
#>     0     1     1
```
