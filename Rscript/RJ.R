# Commands/Functions that might be useful

library(tidyverse)

#-----------------------------------------------------------
# Experiment
  # to input density commands such as dnorm() as an input
experiment <- function(f,x){
  y <- function(x) f(x)
  y(x) # now y is f whatever f is
}

#Checking
experiment(dnorm, 0)
dnorm(0)
experiment(dnorm, 0) == dnorm(0)  #TRUE

ex <- experiment( function(x) dbeta(x, shape1 = 0.2, shape2 = 0.3), 0.4)
dbeta(0.4, shape1 = 0.2, shape2 = 0.3)
ex == dbeta(0.4, shape1 = 0.2, shape2 = 0.3)  #TRUE

#-----------------------------------------------------------
# uniform distribution
  runif()

#-----------------------------------------------------------
# Integrating
  #integrate(f, lower, upper)
  integrate(dnorm, -1.96, 1.96) # approx 0.95
  integrate(dnorm, -Inf, Inf) # approx 1

  # just get the value with $value
  integrate(dnorm, -Inf, Inf)$value  #1
  
#-----------------------------------------------------------  
# Creating g(x) = cf(x)
  # c is the normalizing constant
  g_func <- function(input_func, left_b, right_b, x){
  # docstring:
      # input_func is the underlying density function
      # left_b and right_b are lower and upper bound, respectively
      # x is where you want to evaluate
      # after calculating the normalizing constant c  
      # this function simply creates g as described in the paper
    
    c <- integrate(input_func, left_b, right_b)$value #normalizing constant
    g <- c * input_func(x)
    return(g)
  }

# Checking
  integrate(dnorm, -1.96, 1.96)$value*dnorm(2) #0.05129165
  g_func(dnorm, -1.96, 1.96, 2)                #0.05129165

#-----------------------------------------------------------
# Finding derivative EXPRESSION
deriv(~x^2, "x")                       # .grad[, "x"] gives the expression for derivative
D(expression(x^2), "x")                # just outputs the derivative expression
D(expression(dnorm(x)), "x")           # first derivative for normal density
D(D(expression(dnorm(x)), "x") , "x")  # second derivative for normal density
D(D(expression(x^2), "x") , "x")       # second derivative, 2, which is correct

# for evaluating
eval(D(expression(dnorm(x)), "x"))


one_two_prime <- function(func_x, var, order){
  # docstring
      # this function tries to find the expression for 
      # first derivative or second derivative
  # func_x is the expression that we want to find the derivative of
      # eg expression(func_x), expression(x^2)
  # var will be a string input specifying the variable that we are taking 
  # derivative of
  # order will be either 1 or 2
      # 1 if first derivative
      # 2 if second derivative
  
  if (order < 1 || order > 2) stop("order has to be 1 or 2")
  if (order == 1) D( func_x, var )
  else D( D( func_x, var ), var )
}

# Checking
one_two_prime( expression(dnorm(x)), "x", 0)  # error bc order = 0
one_two_prime( expression(dnorm(x)), "x", 1)  # correct expression
one_two_prime( expression(dnorm(x)), "x", 2)  # correct expression
one_two_prime( expression(dnorm(x)), "x", 3)  # error bc order = 3

one_two_prime( expression(x), "x", 1)         # correct: 1
one_two_prime( expression(x^2), "x", 2)       # correct: 2

#-----------------------------------------------------------
# To EVALUATE the differentiated expressions
x <- 3
one_two_prime( expression(dnorm(x)), "x", 2 ) %>% 
  eval()                         #0.03545479
-(dnorm(3) - 3 * (3 * dnorm(3))) #0.03545479 same as above

# doesn't have to be x and x can be a vector
ha <- c(4, 5, 6)
one_two_prime( expression(dnorm(ha)), "ha", 2 ) %>% eval()


eval_prime <- function(func_x, eval_pts){
  #docstring
    # this function evaluates original function as well as its derivative
    # at specified locations
    # func_x is the function
    # eval_pts is where we want to evaluate
    # first element corresponds to original function
    # second element corresponds to the second function
  
  fprime <- D(func_x, "x") # expression for first derivative
  x <- eval_pts            # points where I want to evaluate
  fprime_eval <- fprime %>% eval() # evaluating first derivative
  f_eval <- func_x %>% eval()      # evaluating original function
  return(list(f_eval, fprime_eval))
}


#Checking
#check with normal distribution
check <- eval_prime(expression(dnorm(x)), c(1,5))
check[[1]]    #2.419707e-01 1.486720e-06     # original function
dnorm(c(1,5)) #2.419707e-01 1.486720e-06

# first derivatives evaluations
  # first derivative expression: -(x * dnorm(x))
check[[2]]                          #-2.419707e-01 -7.433598e-06 
c(-(1 * dnorm(1)), -(5 * dnorm(5))) #-2.419707e-01 -7.433598e-06


#-----------------------------------------------------------
# Log-Concavity Check function
  # by numerically computing the second derivative of the function
  # https://en.wikipedia.org/wiki/Second_derivative <- under Limit section

# Expressions 
rlang::get_expr(body(dnorm))
rlang::get_expr(body(dbeta))

# log wrapper
logfunc <- function(func_x){
  # docstring
    # simply 'log's the function that you want
    # func_x is the function that you want to 'log'
  log_func <- function(x) {}
  old_f <- rlang::get_expr(body(func_x))
  body(log_func) <- rlang::get_expr(quo(log(!!old_f)))
  return(log_func)
}
# Checking
logfunc(function(x) dnorm(x))
beta_dis <- logfunc(function(x) dbeta(x, shape1 = 0.5, shape2 = 0.5))
beta_dis

# Function that numerically calculates 2nd derivative value
numeric_sec_deri <- function(f, x){
  # docstring
    # some functions are hard to derive the second derivative
    # so this calculates second derivative value numerically
    # f is the underlying function
    # x is where we want to evaluate
    # below we are taking limit as h goes to 0
  h <- .Machine$double.eps^(1/4)
  numerator <- f(x+h) - 2*f(x) + f(x-h)
  denominator <- h^2
  val <- numerator / denominator
  return(val)
}
# Checking
numeric_sec_deri(dnorm, x = 2)                           # 0.1619729
numeric_sec_deri(function(x) dnorm(x), x = 2)            # 0.1619729
x <- 2
one_two_prime( expression(dnorm(x)), "x", 2 ) %>% eval() #0.1619729 same as above
numeric_sec_deri(function(x) dnorm(x), x = 1:3) # x can be a vector

# Can find the second derivative 
# even with functions that involve MORE than ONE argument: 
# dbeta() has shape1 and shape2 arguments that need to be provided
numeric_sec_deri(function(x) dbeta(x, shape1 = 0.5, shape2 = 0.5), c(0.2, 0.3))


# Log Concavity checker <- this is the main one
logconcav_check <- function(func, x){
  # docstring
    # This function finds out numerical value of the second derivative
    # at x of the func 
    # if all of the numerical values are negative, then this function
    # tells us that the func is log concavie. vice versa
  result <- numeric_sec_deri(f=logfunc(func), x = x)
  if (max(result) < 0) {     # if max is less than zero, then all are
    cat("The Input Function is Log-Concave")    
  } else {                   # case one or more second derivative evaluation is positive
    cat("The Input Function is NOT Log-Concave.")
  }
  #return(result)
}

#Checking
    # log concave case
    logconcav_check(function(x) dnorm(x), c(-3:18))
    curve(log(dnorm(x)))   # graph looks concave
    # NOT log concave case
    logconcav_check(function(x) dbeta(x, shape1 = 0.2, shape2 = 0.3), c(0.1, 0.2))
    curve(log(dbeta(x, shape1 = 0.2, shape2 = 0.3) )) # definitely not concave
#-----------------------------------------------------------

#Tangent Line from Smoothing Spline Fit
tang_line <- function(x, func_x, x1){
  # docstring:
      # x is the domain we are considering
      # func_x is function of x, ie, y
      # x1 is where we want tangent point to be
      # this function draws a tangent line at x1
  
  y <- function(x) func_x(x) 
  y <- y(x)
  func <- plot(x, y)             # first plot
  spl <- smooth.spline(y~x)      # create spline
  lines(spl, col = "royalblue" ) # draw spline on the plot
  
  x_new <- x1                    # tangent point
  
  y_val <- predict(spl, 
                   x = x_new,    # gives you y-val at x_new 
                   deriv = 0)  
  
  slope <- predict(spl, 
                   x = x_new,    # gives you slope at x_new
                   deriv = 1)    # notice deriv = 1 to get slope
  
  y_int <- y_val$y - (slope$y * x_new) # calculate y intercept of tangent line
  
  x_int <- - y_int / slope$y
  
  # identify on the plot where the tangent point is
  points(y_val, col="yellow", pch=15)
  
  # Draw the tangent line
  lines(x, y_int + slope$y * x, col="red")
  
  return(list(func,        # plot
              y_int,       # y-intercept
              x_int))      # x-intecept
}

#Checking
x <- seq(-2, 2, length.out = 30)
tang_line(x, function(x) dnorm(x), 0.8)
tang_line(x, function(x) x^3, 0.3)

# works with beta distribution 
x <- seq(0.1, 0.99, length.out = 30)
tang_line(x, function(x) dbeta(x, shape1 = 0.2, shape2 = 0.3), 0.3)

#-----------------------------------------------------------
# Segment Joiner
# input should be 
    # function
    # tangent points
    # endpoints of each tangent segments

segment <- function(f, x_vec, end_point){
  # docstring
    # f is the underlying function
    # x_vector contains the tangent points
    # end_point is the endpoints of each tangent line
      # end_point corresponds to Z_j
    # this function connects each tangent line
  
  y <- function(x) f(x)                 # y becomes the underlying function
  y_pt <- y(x_vec)                      # y-coordinates for the x_vec
  func <- plot(x_vec, y_pt)             # first plot the graph
  spl <- smooth.spline(y_pt ~ x_vec)    # create spline
  lines(spl, col = "royalblue", lwd = 3, pch = 3 )   # draw spline on the plot
  
  #calculate slope at each x_vec's element
  slope <- predict(spl, 
                   x = x_vec,    # gives you slope at x_vec
                   deriv = 1) %>% .$y # $y contains the slope elements
  
  # end point calculator
  segment_ender <- function(e){
    # goes through each element of x_vec and calculates y-coordinate of each segment
    
    end <- numeric(length(slope))
    for (i in 1:length(slope)){ #y = slope*(x - x_0) + y_0  <- using simple line equation
      end[i] <- slope[i]*(e[i]-x_vec[i]) + y(x_vec[i]) 
    }
    return(end)
  }
  e <- segment_ender(e=end_point)
  
  # Draw the tangent lines
  lines(end_point,   # end-points of each tangent segment
        e,           # corresponding y coordinates
        lwd = 1.5)
}

# Checking
# with random polynomial
segment(function(x) -x^2+x+1, 
        seq(-1, 1, length.out = 30), 
        end_point = seq(-1.1, 0.99, length.out = 30))
# with sin(x)
segment(function(x) sin(x), 
        seq(-1, 3, length.out = 30), 
        end_point = seq(-1.1, 2.99, length.out = 30))

# with beta distribution
segment(function(x) dbeta(x, shape1 = 0.2, shape2 = 0.3), 
        seq(0.01, 0.99, length.out = 30), 
        end_point = seq(0.01, 0.99, length.out = 30))

#-----------------------------------------------------------
# Function for Z_j
#assume h and hprime are given; they are functions

Z_j <- function(h , hprime, x){
  # docstring
  # Tries to locate z_j's which are points where upper tangent segments
  # intersect with each other
  n <- length(x)
  k <- n - 1
  z <- numeric(n)
  for (i in 1:k){
    numerator <- h(x[i+1]) - h(x[i]) - x[i+1] * hprime(x[i+1]) + x[i] * hprime(x[i])
    denominator <- hprime(x[i]) - hprime(x[i+1])
    z[i] <- numerator/denominator
  }
  return(z)
}

# Checking
Z_j(dnorm, dnorm, c(1,2,3,4,5)) #-0.2872169  0.9105745  1.9688623  2.9887662  0.0000000
Z_j(function(x) dnorm(x), function(x) -(x * dnorm(x)), c(1,2,3,4,5))


# Vectorized way 

Z_j <- function(h , hprime, x){
  # docstring
      # Tries to locate z_j's which are points where upper tangent segments
      # intersect with each other
      # h (=ln(g(x))) is the original function and hprime is the first derivative of h
      # x are the sampled points
  
  h_e <- h(x)                # evaluate original function at x
  h_e_1 <- h_e[-1]           # discard first element 
  
  hprime_e <- hprime(x)      # evaluate first derivative function at x
  hprime_e_1 <- hprime_e[-1] # discard first element 
  
  x_1 <- x[-1]               # discard first element 
  
  n <- length(x)
  z <- numeric(n)            # where result will be stored
  
  numerator <- h_e_1 - h_e - x_1 * hprime_e_1 + x * hprime_e
  denominator <- hprime_e - hprime_e_1
  
  z <- numerator/denominator # formula for z
  return(z)
}

Z_j(dnorm, function(x) -(x * dnorm(x)), c(1,2,3,4,5)) # -0.2872169  0.9105745  1.9688623  2.9887662  0.9999174
#Same answer as the FIRST non-vectorized Z_j except for the last element


check <- eval_prime(expression(dnorm(x)), c(1,5))
#-----------------------------------------------------------
# Function for U_k

U_k <- function(h, hprime, x_vec, x_0){
  # docstring
      # h is the original underlying function
      # hprime is the first derivative of h
      # x_vec contains sampled x
      # x_0 is where we want to evaluate U_k
  h_e <- h(x_vec)
  hprime_e <- hprime(x_vec) 
  h_e + (x_0 - x_vec) * hprime_e
}
# Checking
U_k(dnorm, function(x) -(x * dnorm(x)), c(1,2,3,4,5), 3)

#-----------------------------------------------------------
# Function for S_k

#-----------------------------------------------------------
# Function for l_k
l_k <- function(h, x_vec, x_0){
  # docstring
      # h is the original underlying function
      # x_vec contains the sampled x
      # x_0 is where we want to evaluate l_k
  
  x_vec_1 <- x_vec[-1]
  
  h_e <- h(x_vec)
  h_e_1 <- h(x_vec_1)
  
  num <- (x_vec_1 - x_0) * h_e + (x_0 - x_vec) * h_e_1
  den <- x_vec_1 - x_vec
  
  l <- num / den
  return(l)
}
# Checking 
l_k(dnorm, 1:4, 3)


#-----------------------------------------------------------
# CDF function
cdf <- function(func_x, x){
  # docstring 
    # find the area under func_x and to the left of x
    # func_x is what we are integrating
  f <- function(x) func_x(x)
  v <- integrate(f, -Inf, x)$value
  return(v)
}
cdf(function(x) dnorm(x), 0)



#-----------------------------------------------------------
# Initializing Step

is.infinite(Inf)  # TRUE
is.infinite(-Inf) # TRUE
is.finite(-Inf)   # FALSE

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
#Checking
numeric_first_d(function(x) dnorm(x), 1) # -0.2419707
eval_prime(expression(dnorm(x)), 1)[2]   # -0.2419707  # same as above
numeric_first_d(function(x) dnorm(x), -100) # 0, because graph is essentially flat


starting_x <- function(func_x, min_x, max_x){
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

# Checking
starting_x(function(x) dnorm(x), min_x = 1, max_x=2)
starting_x(function(x) dnorm(x), min_x= -Inf, max_x=2)
starting_x(function(x) dnorm(x), min_x= -1, max_x=Inf)
starting_x(function(x) dnorm(x), min_x= -Inf, max_x=Inf)
starting_x(function(x) dnorm(x), min_x= NA, max_x=2)
starting_x(function(x) dnorm(x), min_x= -1, max_x=NA)
starting_x(function(x) dnorm(x), min_x= NA, max_x=NA)


#-----------------------------------------------------------
# THIS WAS UNSUCCESSUFL 
#Initializing abscissae

#Sampling X_1, ..., X_k

# at x_1, the derivative has to be positive if unbounded on the left
# at x_k, the derivative has to be negative if unbounded on the right

# input function will be given: ex

sorting <- function(x){
  # docstring
    # this function orders element from smallest to the largest AFTER
    # deleting Inf and -Inf
    # input x is a discrete sequence of numbers 
  sort_x <- sort(x) # from smallest to largest
  
  while (sort_x[1] < .Machine$double.xmin){
    sort_x <- sort_x[-1]
  }
  
  while (sort_x[length(sort_x)] > .Machine$double.xmax){
    sort_x <- sort_x[-length(sort_x)]
  }
  return(sort_x)
}

# Checking
sorting(c(-Inf, -Inf, 3, 2, 6, -1, Inf, Inf))


absci <- function(func_x, xmin){
    avg <- mean(xmin, na.rm = T)
    func_x(avg)
    var <- 3
    foo <- one_two_prime(expression(func_x(x)), "x", 1)
    return(foo)
}




eval(one_two_prime(expression(dnorm(x)), "x", 1)) %>% eval(4)
absci(dnorm, 2)
dnorm(3)

func <- dnorm
boo <- 32
one_two_prime(expression(dnorm(boo)), "boo", 1) %>% eval()


one_two_prime(expression(x^2), "x", 1) %>% eval()
