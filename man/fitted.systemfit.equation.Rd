\name{fitted.systemfit.equation}
\alias{fitted.systemfit.equation}
\title{Fitted values}

\description{
   This function extracts the fitted values of one equation
   from an object returned by \code{\link{systemfit}}.
}

\usage{
   \method{fitted}{systemfit.equation}( object, \dots )
}

\arguments{
  \item{object}{an object of type \code{systemfit.equation}.}
  \item{\dots}{other arguments.}
}

\value{
   \code{fitted.systemfit.equation} returns the vector of fitted values.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{
   \code{\link{systemfit}}, \code{\link{fitted.systemfit}}
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

## print the fitted values of the first equation
fitted( fitols$eq[[1]] )

## print the fitted values of the second equation
fitted( fitols$eq[[2]] )
}

\keyword{models}



