test_that("calc_tau", {
  ms = seq(-1, 1, by = 0.01)
  expect_equal(calc_tau(ms), 1)

  ms = seq(-1, 2, by = 0.01)
  expect_lt(calc_tau(ms), 1)

  ms = numeric(100)
  expect_equal(calc_tau(ms), 0)
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
