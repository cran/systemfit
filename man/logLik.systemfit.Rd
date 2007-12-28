\name{logLik.systemfit}
\alias{logLik.systemfit}
\title{Log-Likelihood value of systemfit object}

\description{
   This method calculates the log-likelihood value
   of a fitted object returned by \code{\link{systemfit}}.
}

\usage{
   \method{logLik}{systemfit}( object, ... )
}

\arguments{
   \item{object}{an object of class \code{systemfit}.}
   \item{...}{currently not used.}
}

\value{
   A numeric scalar (the log-likelihood value) with 2 attributes:
   \code{nobs} (total number of observations in all equations) and
   \code{df} (number of free parameters, i.e. coefficients
   + elements of the residual covariance matrix).
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{
   \code{\link{systemfit}}, \code{\link{logLik}}
}


\examples{
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
system <- list( demand = eqDemand, supply = eqSupply )

## perform a SUR estimation
fitsur <- systemfit( system, "SUR", data = Kmenta )

## residuals of all equations
logLik( fitsur )
}

\keyword{models}



