#malvika


min = -Inf

max = Inf

g <- function(x) dnorm(x)

h <- function(x){
  log(g(x))
}

## derivative of h
dh = function(y)
{
  fderiv(h,y,1)
}






z <- function(support)
{
  x0 <- head(support, n=-1)
  x1 <- tail(support, n=-1)
  zed <- x0 + (h(x0) - h(x1) + (x1 - x0)*dh(x1)) / (dh(x1) - dh(x0))
  return(zed)	
}



support <- c(-2,1,2)

z(support)



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


##u_k is totally fine. 


##now the sampling density

y <- 2
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
  
 

plus.cdf(1.5, support)




## sample from the s_k density
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

s_k_sample(support)

  
