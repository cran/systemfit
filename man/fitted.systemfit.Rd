\name{fitted.systemfit}
\alias{fitted.systemfit}
\alias{fitted.systemfit.equation}
\title{Fitted values}

\description{
   These functions extract the fitted values
   from an object returned by \code{\link{systemfit}}.
}

\usage{
   \method{fitted}{systemfit}( object, \dots )

   \method{fitted}{systemfit.equation}( object, \dots )
}

\arguments{
   \item{object}{an object of class \code{systemfit} or
      \code{systemfit.equation}.}
   \item{\dots}{other arguments.}
}

\value{
   \code{fitted.systemfit} returns a data.frame of all fitted values,
   where each column contains the fitted values of one equation.

   \code{fitted.systemfit.equation} returns a vector of the fitted values
      of a single equation.
}

\author{Arne Henningsen \email{arne.henningsen@googlemail.com}}

\seealso{
   \code{\link{systemfit}}, \code{\link{fitted}}
}


\examples{
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
system <- list( demand = eqDemand, supply = eqSupply )

## perform OLS on each of the equations in the system
fitols <- systemfit( system, data = Kmenta )

## all fitted values
fitted( fitols )

## fitted values of the first equation
fitted( fitols$eq[[1]] )

## fitted values of the second equation
fitted( fitols$eq[[2]] )
}

\keyword{models}



