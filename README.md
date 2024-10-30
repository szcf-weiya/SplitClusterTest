# SplitClusterTest

An R package for 

> Wang, L., Lin, Y., & Zhao, H. (2024). False Discovery Rate Control via Data Splitting for Testing-after-Clustering (arXiv:2410.06451). arXiv. <https://doi.org/10.48550/arXiv.2410.06451>

The proposed approach addresses the double-dipping issue in testing-after-clustering tasks, particularly in single-cell data analysis, where the same data is used both for clustering (to identify cell types) and for testing (to select differentially expressed genes), which can inflate false positives.

![dd](https://github.com/user-attachments/assets/e5383503-2e4d-45d0-adff-77f3a0f82899)


## :hammer_and_wrench: Installation

For single-cell data analysis, this package depends on [Seurat](https://github.com/satijalab/Seurat). Since there are potential compatibility issue across different major versions of Seurat. We recommend [Seurat v4](https://github.com/satijalab/seurat/tree/v4.4.0), and it can be installed via the following command:

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
