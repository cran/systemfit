\name{logLik.systemfit}
\alias{logLik.systemfit}
\title{Log-Likelihood value of systemfit object}

\description{
   This method calculates the log-likelihood value
   of a fitted object returned by \code{\link{systemfit}}.
}

\usage{
   \method{logLik}{systemfit}( object, residCovDiag = FALSE, ... )
}

\arguments{
   \item{object}{an object of class \code{systemfit}.}
   \item{residCovDiag}{logical.
      If this argument is set to \code{TRUE},
      the residual covaraince matrix
      that is used for calculating the log-likelihood value
      is assumed to be diagonal,
      i.e. all covariances are set to zero.
      This may be desirable for models
      estimated by OLS, 2SLS, WLS, and W2SLS.}
   \item{...}{currently not used.}
}

\details{
   The residual covariance matrix
   that is used for calculating the log-likelihood value 
   is calculated based on the actually obtained (final) residuals
   (not correcting for degrees of freedom).
   In case of systems of equations with unequal numbers of observations,
   the calculation of the residual covariance matrix
   is only based on the residuals/observations
   that are available in all equations.
}

\value{
   A numeric scalar (the log-likelihood value) with 2 attributes:
   \code{nobs} (total number of observations in all equations) and
   \code{df} (number of free parameters, i.e. coefficients
   + elements of the residual covariance matrix).
}

\author{Arne Henningsen \email{arne.henningsen@googlemail.com}}

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



