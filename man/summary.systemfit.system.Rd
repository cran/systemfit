\name{summary.systemfit.system}
\alias{summary.systemfit.system}

\title{summary.systemfit.system}

\description{
  This function returns a summary of the system of equations.
}

\usage{
summary.systemfit.system(object,...)
}

\arguments{
  \item{object}{an object of type \code{systemfit.system}.}
  \item{...}{not used by user.}
}

\value{
  \code{summary.systemfit.system} returns an object of type \code{systemfit.system}.
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

## perform 3SLS on each of the equations in the system
fit3sls <- threestage.systemfit( system, inst, labels, kmenta )

## print the results
print( fit3sls )



}
\keyword{models}



