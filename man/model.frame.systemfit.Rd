\name{model.frame.systemfit}
\alias{model.frame.systemfit}
\alias{model.frame.systemfit.equation}
\title{Extracting the Data of a systemfit Object}

\description{
   These functions return the data used by \code{\link{systemfit}}
   to estimate a system of equations.
}

\usage{
   \method{model.frame}{systemfit}( formula, \dots )

   \method{model.frame}{systemfit.equation}( formula, \dots )
}

\arguments{
   \item{formula}{an object of class \code{systemfit} or
      \code{systemfit.equation}.}
   \item{\dots}{currently ignored.}
}

\value{
   \code{model.frame.systemfit} returns a simple data frame
      (without a 'terms' attribute) that contains all variables
      used to estimate the entire system of equations.

   \code{model.frame.systemfit.equation} returns a model frame
      (data frame with a 'terms' attribute) that contains
      all variables used to estimate the respective equation.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{
   \code{\link{systemfit}}, \code{\link{model.frame}}, and
   \code{\link{model.matrix.systemfit}}
}


\examples{
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
system <- list( demand = eqDemand, supply = eqSupply )

## perform OLS of the system
fitols <- systemfit( system, data = Kmenta )

## data used to estimate the entire system
model.frame( fitols )

## data used to estimate the first equation
model.frame( fitols$eq[[ 1 ]] )
}

\keyword{models}



