### Auxiliary functions
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


### Initializing Step
starting_x1 <- function(func_x, min_x, max_x){

  # First case
  if (is.finite(min_x) && is.finite(max_x)){ # case when max and min of x are given and FINITE
    opt <- optimize(function(x) func_x(x),
                    c(min_x, max_x),
                    maximum = TRUE)$maximum
    left <- as.vector(quantile(c(min_x, opt, max_x), 0.49))
    right <- as.vector(quantile(c(min_x, opt, max_x), 0.51))
    return(c(left, right))

    # Second case
  } else if ( (is.infinite(min_x) == TRUE || is.na(min_x)==TRUE)  # min is missing or infinte
              && is.finite(max_x) == TRUE ) {
    min_x <- -100
    opt <- optimize(function(x) func_x(x),
                    c(min_x, max_x),
                    maximum = TRUE)$maximum
    left <- as.vector(quantile(c(min_x, opt, max_x), 0.49) )
    right <- as.vector(quantile(c(min_x, opt, max_x), 0.51) )
    return(c(left, right))
  } else if ( (is.infinite(max_x) == TRUE || is.na(max_x)==TRUE) # max is missing or infinite
              && is.finite(min_x)==TRUE ) {
    max_x <- 100
    opt <- optimize(function(x) func_x(x),
                    c(min_x, max_x),
                    maximum = TRUE)$maximum
    left <- as.vector(quantile(c(min_x, opt, max_x), 0.49) )
    right <- as.vector(quantile(c(min_x, opt, max_x), 0.51) )
    return(c(left, right))
  } else {
    max_x <- 100
    min_x <- -100
    opt <- optimize(function(x) func_x(x),
                    c(min_x, max_x),
                    maximum = TRUE)$maximum
    left <- as.vector(quantile(c(min_x, opt, max_x), 0.49) )
    right <- as.vector(quantile(c(min_x, opt, max_x), 0.51) )
    return(c(left, right))
  }
}


#' Adaptive Rejection Sampling - STAT 243 Final Project
#'
#' Conduct adaptive rejection sampling for logarithmically concave input function
#'
#' @import AdapSamp
#' @import rlang
#' @import pracma
#'
#' @param n number of samples
#' @param f input function <- needs to be logarithmically concoave
#' @param min minimum of the domain
#' @param max maximum of the domain
#' @param sp starting points - a vector of 2 points (optional)
#'
#' @return A vector of sampled points, AND a density histogram of the sampled points and the true density curve of the input function
#'
#' @examples
#' a <- ars(300, f = function(x) dnorm(x, 2), -Inf, Inf, c(-2, 3))
#' b <- ars(300, f = function(x) dbeta(x,2,2), 0, 1, c(0.1,0.6))
#'
#' # error case because t-distribution is NOT log-concave
#' ars(20, f = function(x) dt(x,3), -100, 100, c(-2,10))
#'
#' @export
ars <- function(n, f, min, max, sp = NA){


  # Input validation for sample size and domain
  if (n <= 0 || !is.numeric(n))
    stop("Please enter positive integers for sample size.")
  if (min >= max || !is.numeric(min) || !is.numeric(max) )
    stop("Please input valid domain.")

  # create normalized density
  g_func <- function(f, x){
    # docstring:
    # f is the underlying density function
    # left_b and right_b are lower and upper bound, respectively
    # x is where you want to evaluate
    # after calculating the normalizing constant c
    # this function simply creates g as described in the paper

    c <- integrate(f, min, max)$value #normalizing constant
    g <- f(x) / c
    return(g)
  }

  # Create the log of the normalized density
  h <- function(x) {log(g_func(f,x))}
  # Create the function for evaluating the derivative of the h function
  dh <- function(x) {pracma::fderiv(h, x, 1)}
  # Check if the user provides starting points
  if(sum(is.finite(sp)) < length(sp)) {sp = starting_x1(f, min, max)}
  # Check if the starting points are in opposite sides of the density's maximum
  if( prod(dh(sp)) > 0 ) {stop("Starting points must be in opposite sides the maximum")}
  sp <- sort(sp)
  # Initialize sample
  sample <- numeric(n)
  #adjust for the uniform density case
  if (prod(dh(sp) == 0) == 1) { sample <- runif(n, min, max); return(sample)} else {

    # Evaluate the upper hull at a specific point
    u_func = function(x, sp)
    {
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

    # Evaluate the lower hull at a specific point
    l_func <- function(x, sp) {
      # docstring
      # function calculates the lower hull piece corresponding to
      # the data sampled from Sk, set to -inf if x<x1 or x>xk
      # sp is the starting points
      # x is the sampled point
      index <- findInterval(x, sp)
      l <- 0
      if (index == 0 | index == length(sp) ){
        l <- -Inf}
      else {
        l<- ((sp[index + 1] - x) * h(sp[index]) + (x - sp[index]) * h(sp[index + 1])) /
          (sp[index + 1] -sp[index])}

      return(l)
    }

    # Calculate the vector of Z's
    z <- function(support)
    {
      x0 <- head(support, n=-1)
      x1 <- tail(support, n=-1)
      zed <- x0 + (h(x0) - h(x1) + (x1 - x0)*dh(x1)) / (dh(x1) - dh(x0))
      return(zed)
    }

    # Store all the upper hull pieces at one place
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

    # Calculate the probability of intervals
    plus.cdf <- function(sp) {
      #docstring
      #for each starting points, we calculate the value of cdf from zi to zi+1.
      #we store this as a list, each element corresponding to the cdf of that interval.
      ##we normalise the cdf by dividing by its sum.

      zed = c(min,z(sp),max)
      l = length(zed) - 1
      cdf_vals = numeric(l)

      for (i in 1:l) {
        x_j = sp[i]
        cdf_vals[i] = exp(h(x_j))/dh(x_j) * ( exp((zed[i+1] - x_j)*dh(x_j)) - exp((zed[i] - x_j)*dh(x_j)) )}

      normaliser <- sum(cdf_vals)
      return(cdf_vals/normaliser)
    }

    #function to generate xstar from s_k
    s_k_sample <- function(sp) {
      #docstring
      #sampling from cdf: we pick the cdf interval from which to sample from with probability equal to its value
      #then we sample a point from that interval uniformly.
      zed <- c(min,z(sp),max)
      probs <- plus.cdf(sp)

      i <- sample(length(sp), size = 1, prob = probs)
      u <- (runif(1) * probs[i])
      #return(y)
      tmp = (u* dh(sp[i]) / exp(h(sp[i]))) + exp( (zed[i] - sp[i])*dh(sp[i]) )
      y <- ifelse(is.nan(log(tmp)), (tmp-1) - (1/2)*(tmp-1)^2 + (1/3)*(tmp-1)^3, log(tmp))
      y =  (y / dh(sp[i])) + sp[i]
      return(y)

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


    ## begin sampling

    for( i in 1:n){

      accept = 0
      check = 0

      # proceed to sample until one point is accepted by the squeezing/rejection tests
      while(!accept){

        u <- runif(1)
        xstar <- s_k_sample(sp)

        # squeezing and rejection tests
        if (u <= exp(l_func(xstar, sp) - u_func(xstar, sp))) {
          accept = 1
        } else if (u <= exp(h(xstar) - u_func(xstar, sp))) {
          accept = 1
          # include xstar in starting points
          sp <- sort(c(sp, xstar))
          # check the log_concavity of the density for updated starting points
          check <- logconcav_check(sp)
          if (!check) { stop("The function isn't log concave. Please input a log concave function")}
        }
      }

      sample[i] = xstar
    }

    hist(sample, breaks = 200, probability = T)
    curve(f, add = T, col = "yellow", lty = 2, lwd = 2)
    return(sample)
  }
}

