
%    $Id: print.systemfit.equation.Rd,v 1.5 2004/02/23 14:23:46 henningsena Exp $

\name{print.systemfit.equation}
\alias{print.systemfit.equation}
\title{Print output of systemfit estimation}

\description{
  This function prints a summary of one equation of the estimated equation system.
}

\usage{
   \method{print}{systemfit.equation}(x, digits=6,...)
}

\arguments{
  \item{x}{an object of type \code{systemfit.equation}.}
  \item{digits}{number of digits to print.}
  \item{...}{not used by user.}
}

\value{
  \code{print.systemfit.equation} returns nothing.
}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com}}

\seealso{\code{\link{systemfit}}, \code{\link{print.systemfit}} and
  \code{\link{summary.systemfit.equation}}}


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
print( fitols$eq[[1]] )
}

\keyword{models}
