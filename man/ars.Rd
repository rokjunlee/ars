% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ars.R
\name{ars}
\alias{ars}
\title{Adaptive Rejection Sampling - STAT 243 Final Project}
\usage{
ars(n, f, min, max, sp = NA)
}
\arguments{
\item{n}{number of samples}

\item{f}{input function <- needs to be logarithmically concoave}

\item{min}{minimum of the domain}

\item{max}{maximum of the domain}

\item{sp}{starting points - a vector of 2 points (optional)}
}
\value{
A vector of sampled points, AND a density histogram of the sampled points and the true density curve of the input function
}
\description{
Conduct adaptive rejection sampling for logarithmically concave input function
}
\examples{
a <- ars(300, f = function(x) dnorm(x, 2), -Inf, Inf, c(-2, 3))
b <- ars(300, f = function(x) dbeta(x,2,2), 0, 1, c(0.1,0.6))

# error case because t-distribution is NOT log-concave
ars(20, f = function(x) dt(x,3), -100, 100, c(-2,10))

}
