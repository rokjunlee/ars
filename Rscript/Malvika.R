# Malvika

#USEFUL REFERENCE:
> rARS
function (n, formula, min = -Inf, max = Inf, sp) 
{
    sp <- sort(sp)
    if (!is.character(formula)) 
        stop("Unsuitable density function.")
    if (n <= 0) 
        stop("Unsuitable sample size.")
    if (min >= max) 
        stop("Unsuitable domain.")
    p <- function(x) {
        eval(parse(text = formula))
    }
    V <- function(x) {
        -log(p(x))
    }
    x_final <- numeric(n)
    for (j in 1:n) {
        Support <- sp
        if (!identical(Support, sort(Support))) 
            stop("Put the supporting points in ascending order.")
        u = 0
        compareprop = -1
        while (u > compareprop) {
            tangent <- fderiv(V, Support, 1)
            crosspoint = numeric(length(Support) + 1)
            crosspoint[1] = min
            crosspoint[length(crosspoint)] = max
            crossvalue = numeric(length(Support) - 1)
            for (i in 1:(length(Support) - 1)) {
                A = matrix(c(tangent[i], -1, tangent[i + 1], 
                  -1), nrow = 2, byrow = T)
                b = c(tangent[i] * Support[i] - V(Support)[i], 
                  tangent[i + 1] * Support[i + 1] - V(Support)[i + 
                    1])
                solve(A, b)
                crosspoint[i + 1] = solve(A, b)[1]
                crossvalue[i] = solve(A, b)[2]
            }
            IntSum <- numeric(length(Support))
            for (i in 1:length(IntSum)) {
                expfun = function(x) {
                  exp(-tangent[i] * (x - Support[i]) - V(Support)[i])
                }
                IntSum[i] = integrate(expfun, crosspoint[i], 
                  crosspoint[i + 1])[[1]]
            }
            rdm <- runif(1)
            cum = c(0, cumsum(IntSum/sum(IntSum)))
            idx <- which(rdm < cumsum(IntSum/sum(IntSum)))[1]
            x_star <- log((rdm - cum[idx] + exp(tangent[idx] * 
                Support[idx] - V(Support)[idx]) * exp(-tangent[idx] * 
                crosspoint[idx])/sum(IntSum)/(-tangent[idx])) * 
                sum(IntSum) * (-tangent[idx])/exp(tangent[idx] * 
                Support[idx] - V(Support)[idx]))/(-tangent[idx])
            u <- runif(1)
            compareprop <- p(x_star)/exp(-tangent[idx] * (x_star - 
                Support[idx]) - V(Support)[idx])
            Support <- sort(c(Support, x_star))
        }
        x_final[j] = x_star
    }
    x_final
}
<bytecode: 0x155df1ce8>
<environment: namespace:AdapSamp>
