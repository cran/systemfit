\name{print.confint.systemfit}
\alias{print.confint.systemfit}
\title{Print confidence intervals of coefficients}

\description{
  This function prints the confidence intervals of the coefficients
  of the estimated equation system.
}

\usage{
   \method{print}{confint.systemfit}( x, digits=3, \dots )
}

\arguments{
  \item{x}{an object of type \code{confint.systemfit}.}
  \item{digits}{number of digits to print.}
  \item{\dots}{other arguments.}
}

\author{Arne Henningsen \email{arne.henningsen@googlemail.com}}

\seealso{\code{\link{systemfit}}, \code{\link{confint.systemfit}}
  and  \code{\link{confint.systemfit.equation}}}

\examples{
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
system <- list( demand = eqDemand, supply = eqSupply )

## perform OLS on each of the equations in the system
fitols <- systemfit( system, data = Kmenta )

## calculate and print the confidence intervals
## of all coefficients
ci <- confint( fitols )
print( ci, digits=4 )

## calculate and print the confidence intervals
## of the coefficients of the second equation
ci2 <- confint( fitols$eq[[2]] )
print( ci2, digits=4 )
}

\keyword{models}
\keyword{regression}



