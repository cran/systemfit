\name{residuals.systemfit.equation}
\alias{residuals.systemfit.equation}
\title{Residuals of systemfit equation}

\description{
   This function extracts the residuals of one equation
   from an object returned by
   \code{\link{systemfit}}.
}

\usage{
   \method{residuals}{systemfit.equation}( object, \dots)
}

\arguments{
  \item{object}{an object of type \code{systemfit.equation}.}
  \item{\dots}{other arguments.}
}

\value{
   \code{residuals.systemfit.equation} returns a vector of residuals.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{
   \code{\link{systemfit}}, \code{\link{residuals.systemfit}}
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

## print the residuals of the first equation
residuals( fitols$eq[[1]] )

## print the residuals of the second equation
residuals( fitols$eq[[2]] )
}

\keyword{models}



