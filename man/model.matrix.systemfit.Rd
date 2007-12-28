\name{model.matrix.systemfit}
\alias{model.matrix.systemfit}
\alias{model.matrix.systemfit.equation}
\title{Construct Design Matrices for Systems of Equations}

\description{
   These functions create design matrices from objects
   returned by \code{\link{systemfit}}.
}

\usage{
   \method{model.matrix}{systemfit}( object, \dots )

   \method{model.matrix}{systemfit.equation}( object, \dots )
}

\arguments{
   \item{object}{an object of class \code{systemfit} or
      \code{systemfit.equation}.}
   \item{\dots}{currently ignored.}
}

\value{
   \code{model.matrix.systemfit} returns a design matrix to estimate
      the specified system of equations.

   \code{model.matrix.systemfit.equation} returns a design matrix to estimate
      the specified formula of the respective equation.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{
   \code{\link{systemfit}}, \code{\link{model.matrix}}, and
   \code{\link{model.frame.systemfit}}
}


\examples{
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
system <- list( demand = eqDemand, supply = eqSupply )

## perform OLS of the system
fitols <- systemfit( system, data = Kmenta )

## design matrix of the entire system
model.matrix( fitols )

## design matrix of the first equation
model.matrix( fitols$eq[[ 1 ]] )
}

\keyword{models}



