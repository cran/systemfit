\name{vcov.systemfit}
\alias{vcov.systemfit}
\alias{vcov.systemfit.equation}
\title{Variance covariance matrix of coefficients}

\description{
   These functions extract the variance covariance matrix of the
   coefficients from an object returned by
   \code{\link{systemfit}}.
}

\usage{
   \method{vcov}{systemfit}( object, \dots )

   \method{vcov}{systemfit.equation}( object, \dots )
}

\arguments{
   \item{object}{an object of class \code{systemfit} or
      \code{systemfit.equation}.}
   \item{\dots}{other arguments.}
}

\value{
   \code{vcov.systemfit} returns the variance covariance matrix
      of all estimated coefficients.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{
   \code{\link{systemfit}}, \code{\link{vcov}}
}

\examples{
data( "Kmenta" )
demand <- consump ~ price + income
supply <- consump ~ price + farmPrice + trend
labels <- list( "demand", "supply" )
system <- list( demand, supply )

## perform OLS on each of the equations in the system
fitols <- systemfit( "OLS", system, labels, data = Kmenta )

## variance covariance matrix of all coefficients
vcov( fitols )

## variance covariance matrix of the coefficients in the first equation
vcov( fitols$eq[[1]] )

## variance covariance matrix of the coefficients in the second equation
vcov( fitols$eq[[2]] )
}

\keyword{models}



