\name{coef.systemfit.equation}
\alias{coef.systemfit.equation}
\title{Coefficients of systemfit equation}

\description{
   This function extracts the coefficients of a single equation
   from an object returned by \code{\link{systemfit}}.
}

\usage{
   \method{coef}{systemfit.equation}( object, \dots )
}

\arguments{
  \item{object}{an object of type \code{systemfit.equation}.}
  \item{\dots}{other arguments.}
}

\value{
   \code{coef.systemfit.equation} returns a vector of the estimated coefficients
      of a single equation.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{
   \code{\link{systemfit}}, \code{\link{coef.systemfit}}
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

## print the coefficients of the first equation
coef( fitols$eq[[1]] )

## print the coefficients of the second equation
coef( fitols$eq[[2]] )
}

\keyword{models}



