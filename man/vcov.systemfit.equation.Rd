\name{vcov.systemfit.equation}
\alias{vcov.systemfit.equation}
\title{Variance covariance matrix of coefficients}

\description{
   This function extracts the variance covariance matrix of the
   coefficients of one equation from an object returned by
   \code{\link{systemfit}}.
}

\usage{
   \method{vcov}{systemfit.equation}( object, \dots)
}

\arguments{
  \item{object}{an object of type \code{systemfit.equation}.}
  \item{\dots}{other arguments.}
}

\value{
   \code{vcov.systemfit.equation} returns the variance covariance matrix
      of the estimated coefficients.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{
   \code{\link{systemfit}}, \code{\link{vcov.systemfit}}
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

## print the coefficients of the first equation
vcov( fitols$eq[[1]] )

## print the coefficients of the second equation
vcov( fitols$eq[[2]] )
}

\keyword{models}



