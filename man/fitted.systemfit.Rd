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
   \code{fitted.systemfit} returns a vector of all fitted values.

   \code{fitted.systemfit.equation} returns a vector of the fitted values
      of a single equation.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{
   \code{\link{systemfit}}, \code{\link{fitted}}
}


\examples{
data( "Kmenta" )
demand <- consump ~ price + income
supply <- consump ~ price + farmPrice + trend
labels <- list( "demand", "supply" )
system <- list( demand, supply )

## perform OLS on each of the equations in the system
fitols <- systemfit( "OLS", system, labels, data = Kmenta )

## all fitted values
fitted( fitols )

## fitted values of the first equation
fitted( fitols$eq[[1]] )

## fitted values of the second equation
fitted( fitols$eq[[2]] )
}

\keyword{models}



