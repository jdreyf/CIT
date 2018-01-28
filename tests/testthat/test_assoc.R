library(cit)
context("assoc")

set.seed(42)
n <- 10
#x & y negatively correlated
x <- 1:n + rnorm(n); y <- n:1 + rnorm(n); z <- rnorm(n)

test_that("assoc two.sided", {
  expect_equal(assoc(x, y), cor.test(x,y, alternative = "two.sided")$p.value)
  expect_equal(assoc(x, y, z), ppcor::pcor.test(x,y,z=z)$p.value)
})

test_that("assoc one.sided", {
  expect_equal(assoc(x, y, alt="less"), cor.test(x,y, alternative = "less")$p.value)
  expect_lt(assoc(x, y, z, alt="less"), assoc(x, y, z, alt="greater"))
  expect_equal(assoc(x, y, z, alt="less"), ppcor::pcor.test(x,y,z=z)$p.value/2)
  expect_equal(assoc(x, y, z, alt="greater"), 1-ppcor::pcor.test(x,y,z=z)$p.value/2)
})
