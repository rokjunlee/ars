library(testthat)


source("Rscript/RJ.R")


# ------functions-------------------------
##########################################

logfunc <- function(func_x){
  # docstring
  # simply 'log's the function that you want
  # func_x is the function that you want to 'log'
  log_func <- function(x) {}
  old_f <- rlang::get_expr(body(func_x))
  body(log_func) <- rlang::get_expr(quo(log(!!old_f)))
  return(log_func)
}

numeric_first_d <- function(f, x){
  # docstring
  # this function numerically calculates first derivative
  # f is the underlying function
  # x is where we want to evaluate
  # below we are taking limit as h goes to 0
  h <- .Machine$double.eps^(1/4)
  numerator <- f(x+h) - f(x)
  denominator <- h
  val <- numerator / denominator
  return(val)
}

numeric_sec_deri <- function(f, x){
  # docstring
  # some functions are hard to derive the second derivative
  # so this calculates second derivative value numerically
  # f is the underlying function
  # x is where we want to evaluate
  # below we are taking limit as h goes to 0
  h <- .Machine$double.eps^0.25
  numerator <- f(x+h) - 2*f(x) + f(x-h)
  denominator <- h^2
  val <- numerator / denominator
  return(val)
}

g_func <- function(f, x){
  # docstring:
  # f is the underlying density function
  # left_b and right_b are lower and upper bound, respectively
  # x is where you want to evaluate
  # after calculating the normalizing constant c  
  # this function simply creates g as described in the paper
  g <- f(x) / c
  return(g)
}

h <- function(x) {log(g_func(f,x))}
dh <- function(x) {numeric_first_d(h, x)}  

u_func = function(x, sp){
  # docstring
  # function calculates the upper hull piece corresponding to 
  # the data sampled from Sk
  # sp is the starting points
  # x is the sampled point
  z <- numeric(length(sp) + 1)
  z[1] <- min
  z[length(z)] <- max
  z_middle <- sort(z(sp))
  z[2:(length(z)-1)] <- z_middle
  index = findInterval(x, z)
  u <- h(sp[index]) + (x - sp[index])*dh(sp[index])
  return(u)
}


l_func <- function(x, sp) {
  # docstring
  # function calculates the lower hull piece corresponding to 
  # the data sampled from Sk, set to -inf if x<x1 or x>xk
  # sp is the starting points
  # x is the sampled point
  index <- findInterval(x, c(min, sp, max))
  l <- 0
  
  if (x <= sp[(length(sp))] & x >= sp[1]){
    l<- ((sp[index] - x) * h(sp[index - 1]) + (x - sp[index - 1]) * h(sp[index])) / 
      (sp[index] -sp[index - 1])} 
  else{
    l<- -Inf
  }
  return(l)
}


z <- function(support){
  x0 <- head(support, n=-1)
  x1 <- tail(support, n=-1)
  zed <- x0 + (h(x0) - h(x1) + (x1 - x0)*dh(x1)) / (dh(x1) - dh(x0))
  return(zed)	
}


u_k = function(y, support) 
{
  u_plus = rep(0, length(y))
  zed = z(support)
  
  piecewise.idx = findInterval(y, c(min, zed, max))
  npieces = length(zed) + 2
  for(pidx in 1:npieces){
    yp = y[piecewise.idx == pidx]
    xx = h(support[pidx]) + (yp - support[pidx])*dh(support[pidx])
    u_plus[piecewise.idx == pidx] = xx
  }
  return(u_plus)
}


plus.cdf = function(x, sp) {
  zed = z(sp)
  
  zlen = length(zed) ##NUMBER OF CDF POINTS
  s_k_cdf = numeric(length(x))
  normaliser = 0
  for(i in 0:zlen) {
    if(i == 0)
    {
      z_lower = -Inf
    } else {
      z_lower = zed[i]
    }
    
    if(i == zlen)
    {
      z_upper = Inf
    } else {
      z_upper = zed[i+1]
    }
  }
}
    
    s_k_sample = function(support)
    {
      zed = z(support)
      sp = plus.cdf(zed, support)
      zpct = sp$required_cdf
      normaliser = sp$normaliser
      ub = c(0, zpct, 1)
      
      unif.samp = runif(1)
      
      fidx = findInterval(unif.samp, ub)
      num.intervals = length(ub) - 1
      zlow = c(min, zed)
      res = rep(0, length(unif.samp))
      for(i in 1:num.intervals)
      {
        ui = unif.samp[ fidx == i ]
        
        if(length(ui) == 0)
        {
          next
        }
        
        ## Invert the  CDF
        yp = support[i]
        zm = zlow[i]
        tmp = (ui - ub[i]) * dh(yp) * normaliser / exp(h(yp)) + exp( (zm - yp)*dh(yp) )
        tmp = yp + log(tmp) / dh(yp)
        res[ fidx == i ] = tmp
      }
      return(res)
  }

    #check logconcave
logconcav_check <- function(sp){ #x is starting points given
      # docstring
      # This function finds out numerical value of the second derivative
      # at x of the function 
      # if all of the numerical values are negative, then this function
      # tells us that the func is log concave. vice versa
      x_1 <- sp[-1]
      x_0 <- sp[-length(sp)]
      result <- dh(x_1) - dh(x_0)
      result[which(is.nan(result))] = -Inf
      l <- 0
      if (max(result) <= 0) {     # if max is less than zero, then all are
        l <- TRUE
      } else {                   # case one or more second derivative evaluation is positive
        l <- FALSE
      }
      return(l)
}    

#initialize variables in order to run tests
min <- -Inf
max <- Inf
f <- function(x) dnorm(x)
c <- integrate(f, min, max)$value #normalizing constant

##########################################
# ------functions-------------------------


# g_func

test_that("g_func function gives the density of the normalized function", {
      #tests 
      # expect result to be double type
     expect_type(g_func(f, c(-2,2)), 'double')
      
      # first derivative of x^2 is 2x and at x=2, we have 4
     expect_equal(g_func(f, c(-2,2)), dnorm(c(-2,2)))
})

# h

test_that("h function gives the density of the log of g_func", {
  #tests 
  # expect result to be double type
  expect_type(h(c(-2,2)), 'double')
  
  # first derivative of x^2 is 2x and at x=2, we have 4
  expect_equal(h(c(-2,2)), log(dnorm(c(-2,2))))
})

# dh
test_that("dh function gives the derivative of the h function", {
  #tests 
  # expect result to be double type
  expect_type(dh(c(-2,2)), 'double')
  
  # first derivative of x^2 is 2x and at x=2, we have 4
  expect_equal(dh(c(-2,2)), numeric_first_d(logfunc(f), c(-2,2)))
})

# numeric_first_d
test_that("function that numerically evaluates first derivative", {
  
  #tests 
    # expect result to be double type
    expect_type(numeric_first_d(function(x) x, 1), 'double')
  
    # first derivative of x^2 is 2x and at x=2, we have 4
    expect_equal(round(numeric_first_d(function(x) x^2, 2)), 4)
})



# numeric_sec_deri
test_that("function that numerically evaluates second derivative", { 
  
  # tests
    # expect result to be double
    expect_type(numeric_sec_deri(function(x) x^2, -14), 'double')
    
    # second derivative of x^2 is just 2
    expect_equal(numeric_sec_deri(function(x) x^2, -14), 2)
    
    # second derivative of -x^2+2*x is  -2
    expect_equal(numeric_sec_deri(function(x) -x^2+2*x, -14), -2)
})

# z, u_func, l_func

test_that("function z, u_func, l_func gives approproiate values", {
  # expect result to be double type
  expect_type(z(c(-2, 1)), 'double')
  
  # expect result to be double type
  expect_type(u_func(0.5, c(-2,2)), 'double')
  
  # expect result to be double type
  expect_type(l_func(0.5, c(-2,2)), 'double')
})

# logconcav_check
test_that("log concavity check for input function", {
  
  #tests
    f <<- function(x) dnorm(x)
    # dnorm is log-concave
    expect_true( logconcav_check(c(1,2,3)) )  
  
    # dnorm is log-concave regardless of mean and standard deviation
    f<<- function(x) dnorm(x, 2, 5)
    expect_true(logconcav_check(c(1,2,3)))  
    
    # log debeta with inappropriate parameters is not log-concave
    f<<- function(x) dbeta(x, shape1 = 0.2, shape2 = 0.3)
    expect_false( logconcav_check(c(0.1, 0.2)) ) 
    
    # log debeta with appropriate parameters is log-concave
    f<<- function(x) dbeta(x, shape1 = 2, shape2 = 2)
    expect_true( logconcav_check(c(0.1, 0.6)) ) 
    
    # log gamma with inappropriate parameters is not log-concave
    f<<- function(x) dgamma(x, 0.2)
    expect_false( logconcav_check(c(1, 5)) ) 
    
    # log debeta with appropriate parameters is log-concave
    f<<- function(x) dgamma(x, 5)
    expect_true( logconcav_check(c(1, 5)) ) 
    
    #log t distribution is not log-concave
    f<<- function(x) dt(x, 5)
    expect_false( logconcav_check(c(1, 5)) ) 
})




