\name{coef.systemfit}
\alias{coef.systemfit}
\alias{coef.systemfit.equation}
\alias{coef.summary.systemfit}
\alias{coef.summary.systemfit.equation}
\title{Coefficients of systemfit object}

\description{
   These functions extract the coefficients from an object returned by
   \code{\link{systemfit}}.
}

\usage{
   \method{coef}{systemfit}( object, \dots )

   \method{coef}{systemfit.equation}( object, \dots )

   \method{coef}{summary.systemfit}( object, \dots )

   \method{coef}{summary.systemfit.equation}( object, \dots )
}

\arguments{
   \item{object}{an object of class \code{systemfit},
      \code{systemfit.equation}, \code{summary.systemfit}, or
      \code{summary.systemfit.equation}.}
   \item{\dots}{other arguments.}
}

\value{
   \code{coef.systemfit} returns a vector of all estimated coefficients.

   \code{coef.systemfit.equation} returns a vector of the estimated coefficients
      of a single equation.

   \code{coef.summary.systemfit} returns a matrix of all estimated coefficients,
      their standard errors, t-values, and p-values.

   \code{coef.summary.systemfit.equation} returns a matrix of the estimated
      coefficients of a single equation, their standard errors, t-values,
      and p-values.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{
   \code{\link{systemfit}}, \code{\link{coef}}
}


\examples{
data( "Kmenta" )
demand <- consump ~ price + income
supply <- consump ~ price + farmPrice + trend
labels <- list( "demand", "supply" )
system <- list( demand, supply )

## perform OLS on each of the equations in the system
fitols <- systemfit( "OLS", system, labels, data = Kmenta )

## all coefficients
coef( fitols )
coef( summary( fitols ) )

## coefficients of the first equation
coef( fitols$eq[[1]] )
coef( summary( fitols$eq[[1]] ) )

## coefficients of the second equation
coef( fitols$eq[[2]] )
coef( summary( fitols$eq[[2]] ) )
}

\keyword{models}



