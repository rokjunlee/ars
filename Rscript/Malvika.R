#malvika:::

ars <- function(rfunc, n, startingpoints,domain){
      ##check validity of starting points 
      startingpoints <- sort(startingpoints)
      
      #if no starting points provided, can take 
      #x1 = mean - 1s.d
      #x2 = mean + 1s.d
      
      
      ##normalise density
      
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
        
      
      ##checking logconcavity of function
      h <- function(x) {log(rfunc(x))}
      
      ##first derivate test
      r <- domain ##function should be monotonically decreasing
      if (!fderiv(h, fderiv)==sort(h, decreasing = TRUE){print("function isnt log-concave")
      }
      
      ##second derivative test 
      ## test is to check second derivate is non-positive
      if (fderiv(h, fderiv, 2) > 0) {(print("function isnt log concave"))}
          
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
      
   
      
      


  
  
  
      
      
      
      
      
      
