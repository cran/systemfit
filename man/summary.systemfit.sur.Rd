
% 	$Id: summary.systemfit.sur.Rd,v 1.3 2002/12/18 06:40:08 hamannj Exp $	

\name{summary.systemfit.sur}
\alias{summary.systemfit.sur}
\title{summary.systemfit.sur}

\description{
  This function returns a summary of the system of equations.
}

\usage{
summary.systemfit.sur(object,...)
}

\arguments{
  \item{object}{an object of type \code{sur.systemfit}.}
  \item{...}{not used by user.}
}

\value{
  \code{summary.systemfit.sur} returns an object of type \code{systemfit.sur}.
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

## perform SUR on each of the equations in the system
fitsur <- sur.systemfit( system, inst, labels, kmenta )

## print the results
print( fitsur )
}

\keyword{models}



