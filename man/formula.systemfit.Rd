\name{formula.systemfit}
\alias{formula.systemfit}
\alias{formula.systemfit.equation}
\title{Model Formulae of systemfit Objects}

\description{
   This method extracts the model formulae
   from fitted objects returned by \code{\link{systemfit}}.
}

\usage{
   \method{formula}{systemfit}( x, ... )
   \method{formula}{systemfit.equation}( x, ... )
}

\arguments{
   \item{x}{an object of class \code{systemfit}.}
   \item{...}{currently not used.}
}

\value{
   \code{formula.systemfit.equation} returns the formula
   of a single equation of a \code{systemfit} object.
   \code{formula.systemfit.equation} returns a list of formulae:
   one formula object for each equation
   of the \code{systemfit} object.
}

\author{Arne Henningsen \email{arne.henningsen@googlemail.com}}

\seealso{
   \code{\link{systemfit}}, \code{\link{formula}}
}


\examples{
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
system <- list( demand = eqDemand, supply = eqSupply )

## perform a SUR estimation
fitsur <- systemfit( system, "SUR", data = Kmenta )

## formula of the second equation
formula( fitsur$eq[[2]] )

## all formulae of the system
formula( fitsur )
}

\keyword{models}



