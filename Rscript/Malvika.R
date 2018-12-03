#malvika:::
ars <- function(rfunc, n, startingpoints,min,max){
      ##check validity of starting points 
      startingpoints <- sort(startingpoints)
      #sanity check:
      
      fderiv(f, start[1]) > 0
      fderiv(f, start[2]) < 0 
      
      ##calculating z
      zfix = function(support) ###a function of the support
      {
        yfo = head(support, n=-1) ##remove the last element
        yf1 = tail(support, n=-1) ##remove the first element
        zfixed = yf0 + (h(yf0) - h(yf1) + (yf1 - yf0)*derivative(yf1)) / (derivative(yf1) - derivative(yf0))
        return(zfixed)	
      }
      
      ######
      #### unnormalized piecewise-linear upper-bound of the log-density
      
      u_plus = function(y, support) 
      {
        uplus = rep(0, length(y)) ##changes after every iteration
        zfixed = zfix(support) ##get the Z's
        
        
        #using findinterval
        #findInterval(5,c(1,4,9,10)) will give us 2
        
        piecewise.idx = findInterval(y, c(min, zfixed, max))
        npieces = length(zfixed) + 2
        for (pidx in 1:npieces){
          y_piece = y[piecewise.idx == pidx]
          x = logfunction(support[pidx]) + (y_piece - support[pidx])*derivate(support[pidx])
          uplus[piecewise.idx == pidx] = x
        }
        return(res)
      }
      
      ##In particular, the above formulation means that we can precompute G+(zi) which means it is only 
      ##necessary to compute the last, non-zero integral for each y.
      
      uplus.cdf = function(vals, support) 
      {
        # equivalently:  integrate(function(z) exp(hplus(z, support)), lower=-Inf, upper = vals)
        
        zfixed = zfix(support)
        
        zlen = length(support)
        cdf_not_normalised = numeric(length(vals))
        normaliser = 0
        for(z_i in 0:zlen) { ##definition of CDF
          if(z_i == 0)
          {
            zm = -Inf
          } else {
            zm = zfixed[zi]
          }
          
          if(z_i == zlen)
          {
            zp = Inf
          } else {
            zp = zfixed[z_i+1]  ##list of Zs to find cumulative probabilities for
          }
          
          y_piece = support[z_i+1]
          delta = exp(log(y_piece))/derivative(y_piece) * ( exp((zp - y_piece)*derivative(y_piece)) - exp((zm - y_piece)*derivative(y_piece)) )
          
          cidx = zm < vals & vals <= zp
          hidx = vals > zp
          
          cdf_not_normalised[cidx] = cdf_not_normalised[cidx] + exp(h(y_piece))/derivative(y_piece) * ( exp((vals[cidx] - y_piece)*derivative(y_piece)) - exp((zm - y_piece)*derivative(y_piece)) )
          cdf_not_normalised[hidx] = cdf_not_normalised[hidx] + ds
          
          normaliser = normaliser + delta
        }
        
        cdf_upperhull = list( 
          cdf_normalised = cdf_not_normalised / normaliser, 
          normaliser = normaliser
        )
        return(cdf_upperhull)
      }
      
      
      ###sampling from S, arbitrary number of samples
      ##inverting realizations from a Unif(0,1) distribution. Using the previous sum-of-integrals formulation for G+,
      #this requires a search across {G+(z1),⋯G+(zk−1)} and then inverting a single integral.
      
      x_star = function(samp.size = 1, support)
      {
        zfixed = zfix(support)
        gp = uplus.cdf(zfixed, support)
        zcdf = gp$cdf_normalised
        normaliser = gp$normaliser
        upper_bound = c(0, zcdf_not_normalised, 1)
        
        uniform_sample = runif(samp.size)
        
        fidx = findInterval(uniform_sample, upper_bound)
        number_of_intervals = length(upper_bound) - 1
        zlow = c(min, zfixed)
        xstar = rep(NaN, length(uniform_sample)) ##1 for us
        for(i in 1:number_of_intervals)
        {
          ui = uniform_sample[ fidx == i ]
          
          if(length(ui) == 0)
          {
            next
          }
         
          ## Invert the gplus CDF
          y_piece = support[i]
          zm = zlow[i]
          tmp = (ui - upper_bound[i]) * derivative(y_piece) * normaliser / exp(h(y_piece)) + exp( (zm - y_piece)*derivative(y_piece) )
          tmp = y_piece + log(tmp) / derivative(y_piece)
          res[ fidx == i ] = tmp
        }
        return(xstar)
      }
   
      
      


  
  
  
      
      
      
      
      
      
