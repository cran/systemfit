
%    $Id: summary.systemfit.equation.Rd,v 1.5 2004/02/23 14:23:47 henningsena Exp $

\name{summary.systemfit.equation}
\alias{summary.systemfit.equation}
\title{Print output of systemfit estimation}

\description{
  This function returns a summary of one equation of the estimated equation system.
}

\usage{
   \method{summary}{systemfit.equation}( object, ... )
}

\arguments{
  \item{object}{an object of type \code{systemfit.equation}.}
  \item{...}{not used by user.}
}

\value{
  \code{summary.systemfit.equation} returns an object of type \code{systemfit.equation}.
}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com}}

\seealso{\code{\link{systemfit}}, \code{\link{summary.systemfit}} and
  \code{\link{print.systemfit.equation}}}


\examples{
\dontrun{library( systemfit )}

data( kmenta )
attach( kmenta )
demand <- q ~ p + d
supply <- q ~ p + f + a
inst <- ~ d + f + a
labels <- list( "demand", "supply" )
system <- list( demand, supply )

## perform OLS on each of the equations in the system
fitols <- systemfit( "OLS", system, labels )

## print the results
summary.systemfit.equation( fitols$eq[[1]] )
}

\keyword{models}
