\name{predict.systemfit.equation}
\alias{predict.systemfit.equation}

\title{Predictions from Equation System Estimation}

\description{
   Returns the predicted values.
}

\usage{
\method{predict}{systemfit.equation}( object, data, ... )
}

\arguments{
  \item{object}{an object of type \code{systemfit.equation}.}
  \item{data}{data frame in which to predict.}
  \item{...}{additional optional arguments.}
}

\value{
   \code{predict.systemfit.equation} returns a vector
   of the predicted values.
}


\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}


\seealso{
   \code{\link{systemfit}}, \code{\link{predict.systemfit}}
   and \code{\link{predict}}
}

\examples{
\dontrun{library( systemfit )}
data( kmenta )
demand <- q ~ p + d
supply <- q ~ p + f + a
labels <- list( "demand", "supply" )
system <- list( demand, supply )

## OLS estimation
fitols <- systemfit("OLS", system, labels, data=kmenta )

## print the predicted values of the first equation
predict( fitols$eq[[1]] )

## print the predicted values of the second equation
predict( fitols$eq[[2]] )
}

\keyword{models}



