\name{coef.systemfit}
\alias{coef.systemfit}
\title{Coefficients of systemfit object}

\description{
   This function extracts the coefficients from an object returned by
   \code{\link{systemfit}}.
}

\usage{
   \method{coef}{systemfit}( object, \dots)
}

\arguments{
  \item{object}{an object of type \code{systemfit}.}
  \item{\dots}{other arguments.}
}

\value{
   \code{coef.systemfit} returns a vector of all estimated coefficients.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{
   \code{\link{systemfit}},  \code{\link{coef.systemfit.equation}}
   and \code{\link{coef}}
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

## print the coefficients
coef( fitols )
}

\keyword{models}



