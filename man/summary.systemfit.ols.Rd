\name{summary.systemfit.ols}
\alias{summary.systemfit.ols}

\title{summary.systemfit.ols}

\description{
  This function returns a summary of the system of equations.
}

\usage{
summary.systemfit.ols(object,...)
}

\arguments{
  \item{object}{an object of type \code{ols.systemfit}.}
  \item{...}{not used by user.}
}

\value{
  \code{summary.systemfit.ols} returns an object of type \code{systemfit.ols}.
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

## perform OLS on each of the equations in the system
fit1sls <- ols.systemfit( system, inst, labels, kmenta )

## print the results
print( fit1sls )
}

\keyword{models}



