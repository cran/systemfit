
% 	$Id: summary.systemfit.threestage.Rd,v 1.3 2002/12/18 06:40:08 hamannj Exp $	

\name{summary.systemfit.threestage}
\alias{summary.systemfit.threestage}
\title{summary.systemfit.threestage}

\description{
  This function returns a summary of the system of equations.
}

\usage{
summary.systemfit.threestage(object,...)
}

\arguments{
  \item{object}{an object of type \code{threestage.systemfit}.}
  \item{...}{not used by user.}
}

\value{
  \code{summary.systemfit.threestage} returns an object of type \code{systemfit.threestage}.
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

## perform 3SLS on each of the equations in the system
fit3sls <- threestage.systemfit( system, inst, labels, kmenta )

## print the results
print( fit3sls )



}
\keyword{models}



