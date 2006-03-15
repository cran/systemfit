
%    $Id: print.systemfit.Rd 173 2006-02-20 09:10:03Z henningsena $

\name{print.systemfit}
\alias{print.systemfit}
\alias{print.systemfit.equation}
\title{Print results of systemfit estimation}

\description{
   These functions print a few results of the estimated equation system.
}

\usage{
   \method{print}{systemfit}( x, digits=6, \dots )

   \method{print}{systemfit.equation}( x, digits=6, ... )
}

\arguments{
   \item{x}{an object of class \code{systemfit} or \code{systemfit.equation}.}
   \item{digits}{number of digits to print.}
   \item{\dots}{other arguments.}
}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com},\cr
   Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de} }

\seealso{\code{\link{systemfit}}, \code{\link{summary.systemfit}}}


\examples{
data( "Kmenta" )
demand <- consump ~ price + income
supply <- consump ~ price + farmPrice + trend
labels <- list( "demand", "supply" )
system <- list( demand, supply )

## perform OLS on each of the equations in the system
fitols <- systemfit( "OLS", system, labels, data = Kmenta )

## results of the whole system
print( fitols )

## results of the first equation
print( fitols$eq[[1]] )

## results of the second equation
print( fitols$eq[[2]] )
}

\keyword{models}



