\name{print.systemfit.system}
\alias{print.systemfit.system}

\title{print.systemfit.system}

\description{
  This function prints a summary of the system of equations.
}

\usage{
print.systemfit.system(x,digits=6,...)
}

\arguments{
  \item{x}{an object of type \code{threestage.systemfit}.}
  \item{digits}{the number of digits to print.}
  \item{...}{not used by user.}
}

\value{
  \code{print.systemfit.system} returns nothing.
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
fit3sls <- ols.systemfit( system, inst, labels, kmenta )

## print the results
print( fit3sls )

}
\keyword{models}



