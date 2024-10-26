# SplitClusterTest

An R package for 

> Wang, L., Lin, Y., & Zhao, H. (2024). False Discovery Rate Control via Data Splitting for Testing-after-Clustering (arXiv:2410.06451). arXiv. <https://doi.org/10.48550/arXiv.2410.06451>
>

## :hammer_and_wrench: Installation

For single-cell data analysis, this package depends on `Seurat`. Since there are potential compatibility issue across different major versions of Seurat. We recommended Seurat v4, and it can be installed via the following command:

```r
devtools::install_github("satijalab/seurat@v4.4.0")
```

then you can install our package via

```r
devtools::install_github("szcf-weiya/SplitClusterTest")
```

## :arrow_right: See also

- The Julia implementation: <https://github.com/szcf-weiya/SplitClusterTest.jl>
- For the comparison between data splitting and data fission, check <https://github.com/szcf-weiya/fission_vs_splitting>.
