\name{residuals.systemfit}
\alias{residuals.systemfit}
\title{Residuals of systemfit object}

\description{
   This function extracts the residuals
   from an object returned by
   \code{\link{systemfit}}.
}

\usage{
   \method{residuals}{systemfit}( object, \dots)
}

\arguments{
  \item{object}{an object of type \code{systemfit}.}
  \item{\dots}{other arguments.}
}

\value{
   \code{residuals.systemfit} returns a matrix of residuals.
   Each column contains the residuals of one equation.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{
   \code{\link{systemfit}}, \code{\link{residuals.systemfit.equation}}
   and \code{\link{residuals}}
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

## print the residuals
residuals( fitols )

}

\keyword{models}



