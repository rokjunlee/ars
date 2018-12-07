library(testthat)


source("Rscript/RJ.R")



# numeric_sec_deri
test_that("function that numerically evaluates second derivative", { 
  
  # tests
    # expect result to be double
    expect_type(numeric_sec_deri(function(x) x^2, -14), 'double')
    
    # second derivative of x^2 is just 2
    expect_equal(numeric_sec_deri(function(x) x^2, -14), 2)
})


# logconcav_check
test_that("log concavity check for input function", {
  
  #tests
    # log-dnorm is log-concave
    expect_true( logconcav_check(function(x) dnorm(x), c(1,2,3)) )  
    # log debeta is NOT log-concave
    expect_false( logconcav_check( function(x) dbeta(x, shape1 = 0.2, shape2 = 0.3), c(0.1, 0.2)) )
})







