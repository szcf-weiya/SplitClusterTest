# Package index

## All functions

- [`calc_acc()`](https://hohoweiya.xyz/SplitClusterTest/reference/calc_acc.md)
  : Calculate the accuracy (FDR, Power, F1 score) if truth is available

- [`calc_inc_rate()`](https://hohoweiya.xyz/SplitClusterTest/reference/calc_inc_rate.md)
  : Calculate the inclusion rate

- [`calc_tau()`](https://hohoweiya.xyz/SplitClusterTest/reference/calc_tau.md)
  : Calculate the cutoff the mirror statistics

- [`cluster_diff()`](https://hohoweiya.xyz/SplitClusterTest/reference/cluster_diff.md)
  : Calculate the difference across two clusters for each feature

- [`dd()`](https://hohoweiya.xyz/SplitClusterTest/reference/dd.md) :
  Wrapper for the naive double-dipping method

- [`debias_symmetry()`](https://hohoweiya.xyz/SplitClusterTest/reference/debias_symmetry.md)
  : Debias the statistics under the null for symmetry

- [`ds()`](https://hohoweiya.xyz/SplitClusterTest/reference/ds.md) : DS
  procedure

- [`ds(`*`<SingleCellExperiment>`*`)`](https://hohoweiya.xyz/SplitClusterTest/reference/ds.SingleCellExperiment.md)
  : DS procedure for Seurat object

- [`ds(`*`<matrix>`*`)`](https://hohoweiya.xyz/SplitClusterTest/reference/ds.matrix.md)
  : DS procedure for Matrix

- [`est.Sigma()`](https://hohoweiya.xyz/SplitClusterTest/reference/est.Sigma.md)
  : estimate the covariance matrix (assuming there are two clusters)

- [`gen_data_normal()`](https://hohoweiya.xyz/SplitClusterTest/reference/gen_data_normal.md)
  : Generate simulation data with two gaussians

- [`mds()`](https://hohoweiya.xyz/SplitClusterTest/reference/mds.md) :
  MDS procedure

- [`mds(`*`<SingleCellExperiment>`*`)`](https://hohoweiya.xyz/SplitClusterTest/reference/mds.SingleCellExperiment.md)
  : Multiple data splitting

- [`mds(`*`<matrix>`*`)`](https://hohoweiya.xyz/SplitClusterTest/reference/mds.matrix.md)
  : MDS procedure for matrix

- [`mds1()`](https://hohoweiya.xyz/SplitClusterTest/reference/mds1.md) :
  Conduct multiple data splitting

- [`mds1_parallel()`](https://hohoweiya.xyz/SplitClusterTest/reference/mds1_parallel.md)
  : Conduct multiple data splitting in parallel

- [`mds2()`](https://hohoweiya.xyz/SplitClusterTest/reference/mds2.md) :
  Aggregate multiple data splitting results, and return a selection set

- [`mirror_stat()`](https://hohoweiya.xyz/SplitClusterTest/reference/mirror_stat.md)
  : Calculate the mirror statistics

- [`myDiffTTest()`](https://hohoweiya.xyz/SplitClusterTest/reference/myDiffTTest.md)
  : Modified Seurat::DiffTTest by extending the output with statistics
  in addition to p-values

- [`myFindMarkers()`](https://hohoweiya.xyz/SplitClusterTest/reference/myFindMarkers.md)
  : Gene expression markers of identity classes

- [`myGLMDETest()`](https://hohoweiya.xyz/SplitClusterTest/reference/myGLMDETest.md)
  : Modified Seurat::GLMDETest by extending the output with statistics
  in addition to p-values

- [`myWilcoxDETest()`](https://hohoweiya.xyz/SplitClusterTest/reference/myWilcoxDETest.md)
  : Modified Seurat::WilcoxDETest by extending the output with
  statistics in addition to p-values

- [`perform_clustering()`](https://hohoweiya.xyz/SplitClusterTest/reference/perform_clustering.md)
  : Perform Clustering on a SingleCellExperiment object

- [`rankSumTestWithCorrelation()`](https://hohoweiya.xyz/SplitClusterTest/reference/rankSumTestWithCorrelation.md)
  :

  Wilcoxon rank sum test (adapted from
  `limma::rankSumTestWithCorrelation`)

- [`sel_inc_rate()`](https://hohoweiya.xyz/SplitClusterTest/reference/sel_inc_rate.md)
  : Perform selection based on inclusion rate

- [`simdata_1ct`](https://hohoweiya.xyz/SplitClusterTest/reference/simdata_1ct.md)
  :

  Datasets Demo synthetic scRNA-seq data with one cell type based on
  `DuoClustering2018::sce_full_Zhengmix4eq()`

- [`simdata_2ct`](https://hohoweiya.xyz/SplitClusterTest/reference/simdata_2ct.md)
  :

  Demo synthetic scRNA-seq data with two cell types based on
  `DuoClustering2018::sce_full_Zhengmix4eq()`
