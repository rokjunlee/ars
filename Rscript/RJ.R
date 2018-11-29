# Commands/Functions that might be useful

#-----------------------------------------------------------
# uniform distribution
  runif()

#-----------------------------------------------------------
# Integrating
  #integrate(f, lower, upper)
  integrate(dnorm, -1.96, 1.96) # approx 0.95
  integrate(dnorm, -Inf, Inf) # approx 1
  
#-----------------------------------------------------------  
# creating g(x) = cf(x)
  # c is the normalizing constant
  g_func <- function(input_func, left_b, right_b, x){
  # docstring:
      # input_func is the density function
      # left_b and right_b are lower and upper bound, respectively
      # x is where you want to evaluate
      # this function simply creates g as described in the paper
      # after calculating the normalizing constant c
    
    c <- integrate(input_func, left_b, right_b)$value #normalizing constant
    g <- c * input_func(x)
    return(g)
  }

# Check 
  integrate(dnorm, -1.96, 1.96)$value*dnorm(2) #0.05129165
  g_func(dnorm, -1.96, 1.96, 2)                #0.05129165

#-----------------------------------------------------------
# finding derivative
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
#Tangent Line from Smoothing Spline Fit
tang_line <- function(x, func_x, x1){
  # docstring:
      # x is the domain we are considering
      # func_x is function of x, ie, y
      # x1 is where we want tangent point to be
      # this function draws a tangent line at x1
  
  func <- plot(x, func_x)        # first plot
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
              x_int))       # x-intecept
}
x <- seq(0,40)
tang_line(x, y, 32)

#-----------------------------------------------------------
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
#-----------------------------------------------------------
