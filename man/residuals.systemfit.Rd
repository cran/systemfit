\name{residuals.systemfit}
\alias{residuals.systemfit}
\alias{residuals.systemfit.equation}
\title{Residuals of systemfit object}

\description{
   These functions extract the residuals
   from an object returned by
   \code{\link{systemfit}}.
}

\usage{
   \method{residuals}{systemfit}( object, \dots )

   \method{residuals}{systemfit.equation}( object, \dots )
}

\arguments{
   \item{object}{an object of class \code{systemfit} or
      \code{systemfit.equation}.}
   \item{\dots}{other arguments.}
}

\value{
   \code{residuals.systemfit} returns a data.frame of residuals,
   where each column contains the residuals of one equation.

   \code{residuals.systemfit.equation} returns a vector of residuals.
}

\author{Arne Henningsen \email{arne.henningsen@googlemail.com}}

\seealso{
   \code{\link{systemfit}}, \code{\link{residuals}}
}


\examples{
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
system <- list( demand = eqDemand, supply = eqSupply )

## perform OLS on each of the equations in the system
fitols <- systemfit( system, data = Kmenta )

## residuals of all equations
residuals( fitols )

## residuals of the first equation
residuals( fitols$eq[[1]] )

## residuals of the second equation
residuals( fitols$eq[[2]] )
}

\keyword{models}



