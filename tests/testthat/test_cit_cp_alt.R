library(cit)
context("cit_cp_alt")

set.seed(42)
n <- 10
#x -> y -> z
#cor(x,z)<0 & cor(y,z)<0
x <- rep(1:0, each=n/2); y <- x + n:1 + rnorm(n); z <- rnorm(n) - y
C <- rnorm(n)

tmp.t <- cit_cp_alt(L=x,G=y,T=z, n.resampl=10, rseed=0)
tmp.tc <- cit_cp_alt(x,y,z,C=C, n.resampl=10, rseed=0)
tmp.l <- cit_cp_alt(L=x,G=y,T=z, n.resampl=10, rseed=0, alternative="less")
tmp.g <- cit_cp_alt(L=x,G=y,T=z, n.resampl=10, rseed=0, alternative="greater")
tmp.lc <- cit_cp_alt(L=x,G=y,T=z, C=C, n.resampl=10, rseed=0, alternative="less")

test_that("cit_cp_alt two.sided", {
  expect_equal(tmp.t, cit.cp(x,y,z, n.resampl=10, rseed=0))
  expect_equal(tmp.tc, cit.cp(x,y,z, C=C, n.resampl=10, rseed=0))
})

test_that("cit_cp_alt one.sided no C", {
  expect_equal(tmp.t["p_TassocL"]/2, tmp.l["p_TassocL"])
  #don't affect tests with G
  expect_equal(tmp.t["p_TassocGgvnL"], tmp.l["p_TassocGgvnL"])
  expect_equal(tmp.t["p_GassocLgvnT"], tmp.l["p_GassocLgvnT"])
  expect_equal(tmp.t["p_LindTgvnG"], tmp.l["p_LindTgvnG"])
  
  expect_equal(1-tmp.t["p_TassocL"]/2, tmp.g["p_TassocL"])
  #don't affect tests with G
  expect_equal(tmp.t["p_TassocGgvnL"], tmp.g["p_TassocGgvnL"])
  expect_equal(tmp.t["p_GassocLgvnT"], tmp.g["p_GassocLgvnT"])
  expect_equal(tmp.t["p_LindTgvnG"], tmp.g["p_LindTgvnG"])
})

test_that("cit_cp_alt one.sided with C", {
  expect_equal(tmp.tc["p_TassocL"]/2, tmp.lc["p_TassocL"])
  #don't affect tests with G
  expect_equal(tmp.tc["p_TassocGgvnL"], tmp.lc["p_TassocGgvnL"])
  expect_equal(tmp.tc["p_GassocLgvnT"], tmp.lc["p_GassocLgvnT"])
  expect_equal(tmp.tc["p_LindTgvnG"], tmp.lc["p_LindTgvnG"])
})
