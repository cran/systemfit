\name{vcov.systemfit}
\alias{vcov.systemfit}
\title{Variance covariance matrix of coefficients}

\description{
   This function extracts the variance covariance matrix of all
   coefficients from an object returned by
      \code{\link{systemfit}}.
}

\usage{
   \method{vcov}{systemfit}( object, \dots)
}

\arguments{
  \item{object}{an object of type \code{systemfit}.}
  \item{\dots}{other arguments.}
}

\value{
   \code{vcov.systemfit} returns the variance covariance matrix
      of all estimated coefficients.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{
   \code{\link{systemfit}}, \code{\link{vcov.systemfit.equation}}
   and \code{\link{vcov}}
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
vcov( fitols )
}

\keyword{models}



