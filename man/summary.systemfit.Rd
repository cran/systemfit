\name{summary.systemfit}
\alias{summary.systemfit}
\alias{summary.systemfit.equation}
\alias{print.summary.systemfit}
\alias{print.summary.systemfit.equation}
\title{Summary of systemfit estimation}

\description{
   These functions create and print summary results
   of the estimated equation system.
}

\usage{
   \method{summary}{systemfit}( object, ... )

   \method{summary}{systemfit.equation}( object, ... )

   \method{print}{summary.systemfit}( x,
      digits = max( 3, getOption("digits") - 1 ), ... )

   \method{print}{summary.systemfit.equation}( x,
      digits = max( 3, getOption("digits") - 1 ), ... )
}

\arguments{
   \item{object}{an object of class \code{systemfit} or
      \code{systemfit.equation}.}
   \item{x}{an object of class \code{summary.systemfit} or
      \code{summary.systemfit.equation}.}
   \item{digits}{number of digits to print.}
   \item{...}{not used by user.}
}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com},\cr
   Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de} }

\seealso{\code{\link{systemfit}}, \code{\link{print.systemfit}}}


\examples{
data( "Kmenta" )
demand <- consump ~ price + income
supply <- consump ~ price + farmPrice + trend
inst <- ~ income + farmPrice + trend
labels <- list( "demand", "supply" )
system <- list( demand, supply )

## perform OLS on each of the equations in the system
fitols <- systemfit( "OLS", system, labels, data = Kmenta )

## results of the whole system
summary( fitols )

## results of the first equation
summary( fitols$eq[[1]] )

## results of the second equation
summary( fitols$eq[[2]] )
}

\keyword{models}



