\name{print.systemfit.sur}
\alias{print.systemfit.sur}

\title{print.systemfit.sur}

\description{
  This function prints a summary of the system of equations.
}

\usage{
print.systemfit.sur(x,digits=6,...)
}

\arguments{
  \item{x}{an object of type \code{sur.systemfit}.}
  \item{digits}{the number of digits to print.}
  \item{...}{not used by user.}
}

\value{
  \code{print.systemfit.sur} returns nothing.
}

\author{Jeff D. Hamann \email{jeff\_hamann@hamanndonald.com}}

\seealso{\code{\link{ols}},\code{\link{sur}},\code{\link{twostage}} and \code{\link{threestage}}}


\examples{
library( systemfit )

data( kmenta )
attach( kmenta )
demand <- q ~ p + d
supply <- q ~ p + f + a
inst <- ~ d + f + a
labels <- list( "demand", "supply" )
system <- list( demand, supply )

## perform sur on each of the equations in the system
fitsur <- sur.systemfit( system, inst, labels, kmenta )

## print the results
print( fitsur )

}
\keyword{models}



