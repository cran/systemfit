
%    $Id: print.systemfit.Rd,v 1.2 2004/02/23 14:23:46 henningsena Exp $

\name{print.systemfit}
\alias{print.systemfit}
\title{Print output of systemfit estimation}

\description{
  This function prints a summary of the estimated equation system.
}

\usage{
   \method{print}{systemfit}(x, digits=6, \dots)
}

\arguments{
  \item{x}{an object of type \code{systemfit}.}
  \item{digits}{number of digits to print.}
  \item{\dots}{other arguments.}
}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com},\cr
   Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de} }

\seealso{\code{\link{systemfit}}, \code{\link{print.systemfit.equation}}
  and \code{\link{summary.systemfit}}}


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

## print the results
print( fitols )
}

\keyword{models}



