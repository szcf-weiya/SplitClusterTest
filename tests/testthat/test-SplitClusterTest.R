test_that("ds", {
  set.seed(1234)
  data = gen_data_normal(80, 100, delta = 10)
  X = data$X
  res = ds(X, q = 0.1)
  acc = calc_acc(res$sel_set, 1:10)
  expect_lte(acc[["fdr"]], 0.1)
  expect_gte(acc[["power"]], 0.9)

  res2 = mds(X, q = 0.1, M = 10)
  acc2 = calc_acc(res2, 1:10)
  expect_lte(acc2[["fdr"]], 0.1)
  expect_gte(acc2[["power"]], 0.9)
})

test_that("ds_pois", {
  set.seed(1234)
  data = gen_data_pois(n = 100, p = 20, prop = 0.5, delta = 10, rho = 0)
  X = data$X
  res = ds(X, q = 0.2)
  acc = calc_acc(res$sel_set, 1:10)
  expect_lte(acc[["fdr"]], 0.2)
  expect_gte(acc[["power"]], 0.9)
})

test_that("calc_tau", {
  ms = seq(-1, 1, by = 0.01)
  expect_gt(calc_tau(ms), 1)

  ms = seq(-1, 2, by = 0.01)
  expect_lt(calc_tau(ms), 1)

  ms = numeric(100)
  expect_gt(calc_tau(ms), 0)
})

test_that("mds", {
  mss = list()
  mss[[1]] = seq(-1, 2, by = 0.01)
  mss[[2]] = seq(-1, 2, by = 0.01)
  names(mss[[1]]) = 1:length(mss[[1]])
  names(mss[[2]]) = 1:length(mss[[2]])
  sel = mds2(mss)
  val.min = mss[[1]][sel][1]
  expect_lt(val.min, 1)
  expect_lt( (length(sel) - 100) / length(sel), 0.05 )

  sel = mds2(mss, tied.method = 2)
  val.min = mss[[1]][sel][1]
  expect_lt(val.min, 1)
  expect_lt( (length(sel) - 100) / length(sel), 0.05 )

  sel = mds2(mss, tied.method = 3)
  val.min = mss[[1]][sel][1]
  expect_lt(val.min, 1)
  expect_lt( (length(sel) - 100) / length(sel), 0.05 )
})

test_that("clustering", {
  tstat = perform_clustering(simdata_2ct$simu_sce)
  nde = sum(abs(tstat) > 2)
  expect_gt(nde, 0)
  expect_lt(nde, 198)
})

test_that("calc_acc", {
  acc = calc_acc(c(), c(1, 2))
  expect_equal(acc[["fdr"]], 0)
  expect_equal(acc[["power"]], 0)
  expect_equal(acc[["f1"]], 0)

  acc = calc_acc(c(1), c(1, 2))
  expect_equal(acc[["fdr"]], 0)
  expect_equal(acc[["power"]], 0.5)
  expect_equal(acc[["f1"]], 2/3)

  acc = calc_acc(c(), c())
  expect_equal(acc[["fdr"]], 0)
  expect_equal(acc[["power"]], 1)
  expect_equal(acc[["f1"]], 1)

  acc = calc_acc(c(1), c())
  expect_equal(acc[["fdr"]], 1)
  expect_equal(acc[["power"]], 0)
  expect_equal(acc[["f1"]], 0)
})

test_that("est.Sigma", {
  x = rbind(matrix(rnorm(100), 50, 2), matrix(rnorm(100) + 4, 50, 2))
  S = est.Sigma(x)
  S0 = diag(2)
  S1 = cov(x)
  expect_lt(norm(S - S0), norm(S1 - S0))
})
