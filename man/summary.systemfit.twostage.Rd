\name{summary.systemfit.twostage}
\alias{summary.systemfit.twostage}

\title{summary.systemfit.twostage}

\description{
  This function returns a summary of the system of equations.
}

\usage{
summary.systemfit.twostage(object,...)
}

\arguments{
  \item{object}{an object of type \code{twostage.systemfit}.}
  \item{...}{not used by user.}
}

\value{
  \code{summary.systemfit.twostage} returns an object of type \code{systemfit.ols}.
}

\author{Jeff D. Hamann \email{jeff\_hamann@hamanndonald.com}}

\seealso{\code{\link{ols}},\code{\link{twostage}} and \code{\link{threestage}}}


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
fit2sls <- twostage.systemfit( system, inst, labels, kmenta )

## print the results
print( fit2sls )



}
\keyword{models}



