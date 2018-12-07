library(testthat)


source("Rscript/RJ.R")


test_that("standardize works with normal input", {
  x <- c(1, 2, 3)
  z <- (x - mean(x)) / sd(x)
  
  expect_equal(standardize(x), z)
  expect_length(standardize(x), length(x))
  expect_type(standardize(x), 'double')
})


test_that("evaluate second derivative", { 
  x <- seq(-2, 3)
  
  expect_equal(numeric_sec_deri(function(x) x^2, 2), 2)
  expect_equal(numeric_sec_deri(function(x) dnorm(x), 2), 2)
})

