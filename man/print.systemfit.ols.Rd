
% 	$Id: print.systemfit.ols.Rd,v 1.3 2002/12/18 06:40:07 hamannj Exp $	

\name{print.systemfit.ols}
\alias{print.systemfit.ols}
\title{print.systemfit.ols}

\description{
  This function prints a summary of the system of equations.
}

\usage{
print.systemfit.ols(x,digits=6,...)
}

\arguments{
  \item{x}{an object of type \code{ols.systemfit}.}
  \item{digits}{number of digits to print.}
  \item{...}{not used by user.}
}

\value{
  \code{print.systemfit.ols} returns nothing.
}

\author{Jeff D. Hamann \email{jeff\_hamann@hamanndonald.com}}

\seealso{\code{\link{ols}},\code{\link{twostage}},\code{\link{sur}}and \code{\link{threestage}}}



\examples{
library( systemfit )

data( kmenta )
attach( kmenta )
demand <- q ~ p + d
supply <- q ~ p + f + a
inst <- ~ d + f + a
labels <- list( "demand", "supply" )
system <- list( demand, supply )

## perform OLS on each of the equations in the system
fit1sls <- ols.systemfit( system, inst, labels, kmenta )

## print the results
print( fit1sls )



}
\keyword{models}



