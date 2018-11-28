# Functions that we need 

# uniform distribution
  runif()

# Integrating
  #integrate(f, lower, upper)
  integrate(dnorm, -1.96, 1.96) # approx 0.95
  
# creating g(x) = cf(x)
  # c is the normalizing constant
  g_func <- function(input_func, left_b, right_b, x){
    c <- integrate(input_func, left_b, right_b)$value #normalizing constant
    g <- c * input_func(x)
    return(g)
  }
  integrate(dnorm, -1.96, 1.96)$value*dnorm(2) #0.05129165
  g_func(dnorm, -1.96, 1.96, 2)                #0.05129165

  
  
