\name{confint.systemfit}
\alias{confint.systemfit}
\alias{confint.systemfit.equation}
\title{Confidence intervals of coefficients}

\description{
   These functions calculate the confidence intervals of the
   coefficients from an object returned by \code{\link{systemfit}}.
}

\usage{
   \method{confint}{systemfit}( object, parm = NULL, level = 0.95,
      useDfSys = NULL, \dots )

   \method{confint}{systemfit.equation}( object, parm, level = 0.95,
      useDfSys = NULL, \dots )
}

\arguments{
   \item{object}{an object of class \code{systemfit} or \code{systemfit.equation}.}
   \item{parm}{not used yet.}
   \item{level}{confidence level.}
   \item{useDfSys}{logical. Use the degrees of freedom of the whole system
      (in place of the degrees of freedom of the single equation)
      to calculate the confidence intervals of the coefficients.
      If it not specified (\code{NULL}), it is set to \code{TRUE}
      if there are parameter restrictions and \code{FALSE} otherwise.}
   \item{\dots}{other arguments.}
}

\value{
   An object of class \code{confint.systemfit}, which is a matrix
   with columns giving lower and upper confidence limits
   for each estimated coefficient. These will be labelled as
   (1-level)/2 and 1 - (1-level)/2 in \% (by default 2.5\% and 97.5\%).
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{
   \code{\link{systemfit}},
   \code{\link{print.confint.systemfit}}, \code{\link{confint}}
}


\examples{
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
system <- list( demand = eqDemand, supply = eqSupply )

## perform OLS on each of the equations in the system
fitols <- systemfit( system, data = Kmenta )

## confidence intervals of all coefficients
confint( fitols )

## confidence intervals of the coefficients of the first equation
confint( fitols$eq[[1]] )

## confidence intervals of the coefficients of the second equation
confint( fitols$eq[[2]] )
}

\keyword{models}
\keyword{regression}



