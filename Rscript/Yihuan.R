#Yihuan

## rARS function with comments 

# Arguments
# n ： Desired sample size;
# formula ： Kernal of the target density;
# min, max ： Domain including positive and negative infinity of the target distribution; 
# sp ：  Supporting set.

# Example 1: Standard normal distribution
# x1 <- rARS(100,"exp(-x^2/2)",-Inf,Inf,c(-2,2))


rARS <- function (n, formula, min = -Inf, max = Inf, sp) 
{
  
  sp <- sort(sp)   ###starting points(?)
  
  ### input validations
  if (!is.character(formula)) 
    stop("Unsuitable density function.")
  if (n <= 0) 
    stop("Unsuitable sample size.")
  if (min >= max) 
    stop("Unsuitable domain.")
  
  ###convert formula from text to actual density (kernal)
  p <- function(x) {
    eval(parse(text = formula))
  }
  
  ###evaluate neg-log of the density( -h(x) )
  V <- function(x) {
    -log(p(x))
  }
  
  x_final <- numeric(n)  #n 0s
  
  
  for (j in 1:n) {
    
    #check support ascending order
    Support <- sp                                  
    if (!identical(Support, sort(Support))) 
      stop("Put the supporting points in ascending order.")
    u = 0
    compareprop = -1
    while (u > compareprop) {
      
      ###first derivative of -h(x) at sp, -h'(x)
      tangent <- fderiv(V, Support, 1)           
      
      ###initialize z
      crosspoint = numeric(length(Support) + 1)   # initialize
      crosspoint[1] = min                         # -inf
      crosspoint[length(crosspoint)] = max        # inf
      crossvalue = numeric(length(Support) - 1)   # initialize
      
      ###calculate z
      for (i in 1:(length(Support) - 1)) {
        A = matrix(c(tangent[i], -1, tangent[i + 1],             ##A = -h'(x1) -1
                                                                 ##     -h'(x2) -1
                     -1), nrow = 2, byrow = T)
        b = c(tangent[i] * Support[i] - V(Support)[i],               ##b = (-h'(x1)*x1 -(-h(x1)),-h'(x2)*x2 -(-h(x2)))
              tangent[i + 1] * Support[i + 1] - V(Support)[i + 1]) 
        solve(A, b)       
        crosspoint[i + 1] = solve(A, b)[1]            ##Z1
        crossvalue[i] = solve(A, b)[2]     
      }
      # crosspoint
      #[1] -Inf    0  Inf
      # crossvalue
      #[1] -2
      
      ### calculate Uk(x) and integral
      IntSum <- numeric(length(Support))   #0 0
      for (i in 1:length(IntSum)) {
        expfun = function(x) {
          exp(-tangent[i] * (x - Support[i]) - V(Support)[i])
        }
        IntSum[i] = integrate(expfun, crosspoint[i], 
                              crosspoint[i + 1])[[1]]     #intsum = 3.694528 3.694528
      }
      
      ###sample x_star 
      rdm <- runif(1)
      cum = c(0, cumsum(IntSum/sum(IntSum)))    #0.0 0.5 1.0 (cdf)
      idx <- which(rdm < cumsum(IntSum/sum(IntSum)))[1]
      x_star <- log((rdm - cum[idx] + exp(tangent[idx] * 
          Support[idx] - V(Support)[idx]) * exp(-tangent[idx] * 
          crosspoint[idx])/sum(IntSum)/(-tangent[idx])) * 
          sum(IntSum) * (-tangent[idx])/exp(tangent[idx] * 
          Support[idx] - V(Support)[idx]))/(-tangent[idx])    
      u <- runif(1)
      compareprop <- p(x_star)/exp(-tangent[idx] * (x_star - 
                                                      Support[idx]) - V(Support)[idx])
      Support <- sort(c(Support, x_star)) #include xstar in Tk
    }
    x_final[j] = x_star #update sample
  }
  x_final
}

#sp is vector of starting points

h <- function(x) log(g(x))
t <- fderiv(h, sp, 1)

u <- h(sp[i]) + (x - sp[i])*t[i]
s <- exp(u) / integrate(exp(u), z0, zk)
l <- ((sp[i+1] - x) * h(sp[i]) + (x - sp[i]) * h(sp[i+1])) / (sp[i+1] - sp[i])
z <- (h(sp[i+1]) - h(sp[i]) - sp[i+1]*t[i+1] + sp[i]*t[i]) / (t[i] - t[i+1])







### basic function structure
ars <- function(n, f, min, max, sp){
  for( i in 1:n){
    
    h <- function(x) log(g(x))
    sp <- sort(sp)
    accept = 0
    sample <- numeric(n)
    
    while(!accept){
      
      hprime <- function(h, sp) {fderiv(h, sp, 1)  #h'(x)}
      
      #find vector of z
        Z_j <- function(h , hprime, sp, x){
          # docstring
          # Tries to locate z_j's which are points where upper tangent segments
          # intersect with each other
          # h (=ln(g(x))) is the original function and hprime is the first derivative of h
          # x are the sampled points
          
          z <- numeric(length(sp) + 1)
          z[1] <- min
          z[length(z)] <- max
          
          h_e <- h(sp)              # evaluate original function at x
          h_e_1 <- h_e[-1]           # discard first element
          h_e_k <- h_e[-length(sp)]          # discard last element 
          
          hprime_e <- hprime(sp)      # evaluate first derivative function at x
          hprime_e_1 <- hprime_e[-1] # discard first element 
          hprime_e_k <- hprime_e[-length(sp)] # discard last element 
          
          sp_1 <- sp[-1]               # discard first element 
          sp_k <- sp[-length(sp)]               # discard last element
          
          numerator <- h_e_1 - h_e_k - sp_1 * hprime_e_1 + sp_k * hprime_e_k
          denominator <- hprime_e_k - hprime_e_1
          
          z <- numerator/denominator    # formula for z
          z <- append(min, z, max)
          
          return(z)
        }
      
      for(k in 1:(length(z) - 2)){
        z[k+1] <- (h(sp[k+1]) - h(sp[k]) - sp[k+1]*t[k+1] + sp[k]*t[k]) / (t[k] - t[k+1])
      }
      
      # check where of z that x belongs to 
      # now assume x belongs to [zj-1,zj]
      u_func <- function(x, j) {
        
        h(sp[j]) + (x - sp[j])*t[j] 
      }
      
      # check where of starting points that x belongs to 
      # now assume x belongs to [xj,xj+1]
      # how to adjust so that ifx<x1 orx>xk,lk(x) is set equal to -inf?
      l_func <- function(x, j) ((sp[j+1] - x) * h(sp[j]) + (x - sp[j]) * h(sp[j+1])) / (sp[j+1] - sp[j])
      
      s <- exp(u_func) / integrate(exp(u_func), z[1], z[length(z)])
      
      
      u <- runif(1)
      #assume sampled xstar from s(how?)
      
      #squeezing and rejection tests
      ifelse(u <= exp(l(xstar) - u(xstar)), accept = 1, ifelse(u <= exp(h(xstar) - u(xstar)),accept = 1))
      #include xstar in sp
      sp <- sort(c(sp, xstar))
    }
    sample[i] = xstar
  }
  return(sample)
  }
  
  
  
  
  
  
  
  
  ###test
  #check logconcave
  
  logconcav_check <- function(func, x){ #x is domain given
    # docstring
    # This function finds out numerical value of the second derivative
    # at x of the func 
    # if all of the numerical values are negative, then this function
    # tells us that the func is log concavie. vice versa
    result <- numeric_sec_deri(f=logfunc(func), x = x)
    l <- 0
    if (max(result) < 0) {     # if max is less than zero, then all are
      l <- TRUE
    } else {                   # case one or more second derivative evaluation is positive
      l <- FALSE
    }
    return(l)
  }
  
  
  # D denotes the domain of f(x)
  # Choose k starting abscissae
  # From pg 342: "In general we have found two starting abscissae (k = 2) to be necessary and sufficient for computational efficiency."
  
  
  
  # need some ways to specify x_1 and x_k
  # better to find x_1 as small as possible as long as derivative is positive and density is positive
  starting_x <- function(func_x, min_x = -Inf, max_x = Inf){
    # docstring
    # divide into 4 cases
    # 1: both min and max are given
    # 2: when max is finite but min is missing or infinite
    # 3: when min is finite but max is missing or infinite
    # 4: otherwise
    
    # func_x is the underlying function
    # min_x and max_x are min and max of f's domain
    
    
    # First case
    if (is.finite(min_x) && is.finite(max_x)){ # case when max and min of x are given and FINITE
      
      return(c(min_x, max_x))  # simplest case: just return min and max values
      
      # Second case  
    } else if ( (is.infinite(min_x) == TRUE || is.na(min_x)==TRUE)  # min is missing or infinte
                && is.finite(max_x) == TRUE ) { # max is finite
      
      if (numeric_first_d(function(x) func_x(x), -.Machine$integer.max) > 0){
        min_x <- -.Machine$integer.max # if derivative is positive then just use -.Machine$integer.max
      } else { # otherwise slowly increase min_x until it becomes positive
        min_x <- -1000  
        while (numeric_first_d(function(x) func_x(x), min_x) <= 0) { 
          print(numeric_first_d(function(x) func_x(x), min_x))
          min_x <- min_x + 5
        }
      }
      return(c(min_x, max_x))
      
      
      # Third case
    } else if ( (is.infinite(max_x) == TRUE || is.na(max_x)==TRUE) # max is missing or infinite
                && is.finite(min_x)==TRUE ) {  # min is finite
      
      if (numeric_first_d(function(x) func_x(x), .Machine$integer.max) < 0){
        max_x <- .Machine$integer.max # if derivative is negative then just use .Machine$integer.max
      } else {  # otherwise slowly decrease max_x until it becomes negative
        max_x <- 1000
        while (numeric_first_d(function(x) func_x(x), max_x) >= 0) {
          print(numeric_first_d(function(x) func_x(x), max_x))
          max_x <- max_x - 5
        }
      }
      return(c(min_x, max_x))
      
      # Fourth case
    } else {  # both min and max are infinite or both missing
      
      if (numeric_first_d(function(x) func_x(x), -.Machine$integer.max) > 0){
        min_x <- -.Machine$integer.max
      } else {
        min_x <- -1000
        while (numeric_first_d(function(x) func_x(x), min_x) <= 0) {
          print(numeric_first_d(function(x) func_x(x), min_x))
          min_x <- min_x + 5
        }
      }
      
      if (numeric_first_d(function(x) func_x(x), .Machine$integer.max) < 0){
        max_x <- .Machine$integer.max
      } else {
        max_x <- 1000
        while (numeric_first_d(function(x) func_x(x), max_x) >= 0) {
          print(numeric_first_d(function(x) func_x(x), max_x))
          max_x <- max_x - 5
        }
      }
      return(c(min_x, max_x))   
      
    }
  }
  
  #calculate the functions u_k(X), s_k(x) and l_k(x)
  # These are just formula <- just function
 

  ### basic function structure
  ars <- function(n, f, min, max, sp = NaN){
    
    if (n <= 0 || !is.numeric(n)) 
      stop("Please enter positive integers for sample size.")
    if (min >= max || !is.numeric(min) || !is.numeric(max) ) 
      stop("Please input valid domain.")
    check <- logconcav_check(f, seq(min, max, by = 0.1))
    if (!check) { print("The function isn't log concave. Please input a log concave function")} else{
      if(sum(is.finite(sp)) < 2) {sp = starting_x(f, min, max)}
      sp <- sort(sp)
      sample <- numeric(n)
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
      
      h <- function(x) {log(g_func(f,x))}
      dh <- function(x) {fderiv(h, x, 1)} 
      
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
      
      
      l_func <- function(x, sp) {
        # docstring
        # function calculates the lower hull piece corresponding to 
        # the data sampled from Sk, set to -inf if x<x1 or x>xk
        # sp is the starting points
        # x is the sampled point
        index <- findInterval(x, c(min, sp, max))
        l <- 0
        
        if (x <= sp[(length(sp))] & x >= sp[1]){
          l<- ((sp[index] - x) * h(sp[index - 1]) + (x - sp[index - 1]) * h(sp[index])) / (sp[index] - sp[index - 1])} else{
            l<- -Inf
          }
        return(l)
        
      }
      
      z <- function(support)
      {
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
      
      plus.cdf <- function(y, support) 
      {
        zed <- c(min, z(support), max)
        p <- findInterval(y, zed)
        l <- length(zed)
        cdf_not_normalised <- numeric(l-1)
        
        for (i in 1:(l-1)){
          cdf_not_normalised[i] <- integrate(function(z) exp(u_k(z, support)),zed[i], zed[i+1])$value
        }
        normaliser <- sum(cdf_not_normalised)
        required_cdf <- (sum(cdf_not_normalised[1:(p-1)]) + integrate(function(z) exp(u_k(z, support)),zed[p],y)$value)  
        
        l <- list(required_cdf = required_cdf/normaliser, normaliser = normaliser)
        return(l)
      }
      
      s_k_sample = function(support)
      {
        zed = z(support)
        sp = sapply(zed, function(x) plus.cdf(x, support), simplify = TRUE)
        zpct = sp[1,]
        norm.const = sp[2][[1]]
        ub = unlist(c(0, zpct, 1))
        
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
          tmp = (ui - ub[i]) * dh(yp) * norm.const / exp(h(yp)) + exp( (zm - yp)*dh(yp) )
          tmp = yp + log(tmp) / dh(yp)
          res[ fidx == i ] = tmp
        }
        return(res)
      }
      
      
      for( i in 1:n){
        
        accept = 0
        
        while(!accept){
          
          u <- runif(1)
          xstar <- s_k_sample(sp)
          
          #squeezing and rejection tests
          if (u <= exp(l_func(xstar, sp) - u_func(xstar, sp))) { 
            accept = 1 
          } else if (u <= exp(h(xstar) - u_func(xstar, sp))) { 
            accept = 1 
            #include xstar in sp
            sp <- sort(c(sp, xstar))}
        }
        sample[i] = xstar
      }
      return(sample)
    }
  }
  ars(10, f = function(x) dnorm(x), -10, 10, c(-2,2))

  
  