\name{terms.systemfit}
\alias{terms.systemfit}
\alias{terms.systemfit.equation}
\title{Model Terms of systemfit Objects}

\description{
   This method extracts the model terms
   from fitted objects returned by \code{\link{systemfit}}.
}

\usage{
   \method{terms}{systemfit}( x, ... )
   \method{terms}{systemfit.equation}( x, ... )
}

\arguments{
   \item{x}{an object of class \code{systemfit}.}
   \item{...}{currently not used.}
}

\value{
   \code{terms.systemfit.equation} returns the model terms
   of a single equation of a \code{systemfit} object.
   \code{terms.systemfit.equation} returns a list of model terms:
   one model term object for each equation
   of the \code{systemfit} object.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{
   \code{\link{systemfit}}, \code{\link{terms}}
}


\examples{
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
system <- list( demand = eqDemand, supply = eqSupply )

## perform a SUR estimation
fitsur <- systemfit( system, "SUR", data = Kmenta )

## model terms of the second equation
terms( fitsur$eq[[ 2 ]] )

## all model terms of the system
terms( fitsur )
}

\keyword{models}

