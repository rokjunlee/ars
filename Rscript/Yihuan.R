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
