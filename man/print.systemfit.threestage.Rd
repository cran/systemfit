
% 	$Id: print.systemfit.threestage.Rd,v 1.3 2002/12/18 06:40:08 hamannj Exp $	

\name{print.systemfit.threestage}
\alias{print.systemfit.threestage}
\title{print.systemfit.threestage}

\description{
  This function prints a summary of the system of equations.
}

\usage{
print.systemfit.threestage(x,digits=6,...)
}

\arguments{
  \item{x}{an object of type \code{threestage.systemfit}.}
  \item{digits}{the number of digits to print.}
  \item{...}{not used by user.}
}

\value{
  \code{print.systemfit.threestage} returns nothing.
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

## perform 2SLS on each of the equations in the system
fit3sls <- ols.systemfit( system, inst, labels, kmenta )

## print the results
print( fit3sls )

}
\keyword{models}



