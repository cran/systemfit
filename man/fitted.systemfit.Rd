\name{fitted.systemfit}
\alias{fitted.systemfit}
\title{Fitted values}

\description{
   This function extracts the fitted values of all equations
   from an object returned by \code{\link{systemfit}}.
}

\usage{
   \method{fitted}{systemfit}( object, \dots )
}

\arguments{
  \item{object}{an object of type \code{systemfit}.}
  \item{\dots}{other arguments.}
}

\value{
   \code{fitted.systemfit} returns the vector of fitted values.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{
   \code{\link{systemfit}}, \code{\link{fitted.systemfit.equation}}
   and \code{\link{fitted}}
}


\examples{
\dontrun{library( systemfit )}

data( kmenta )
attach( kmenta )
demand <- q ~ p + d
supply <- q ~ p + f + a
labels <- list( "demand", "supply" )
system <- list( demand, supply )

## perform OLS on each of the equations in the system
fitols <- systemfit( "OLS", system, labels )

## print the fitted values
fitted( fitols )
}

\keyword{models}



