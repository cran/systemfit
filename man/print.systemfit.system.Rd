
% 	$Id: print.systemfit.system.Rd,v 1.6 2003/11/10 20:19:05 hamannj Exp $	

\name{print.systemfit.system}
\alias{print.systemfit.system}
\title{print.systemfit.system}

\description{
  This function prints a summary of the estimated equation system.
}

\usage{
print.systemfit.system(x, digits=6,...)
}

\arguments{
  \item{x}{an object of type \code{systemfit.system}.}
  \item{digits}{number of digits to print.}
  \item{...}{not used by user.}
}

\value{
  \code{print.systemfit.system} returns nothing.
}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com}}

\seealso{\code{\link{systemfit}}, \code{\link{print.systemfit.equation}}
  and \code{\link{summary.systemfit.system}}}


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



