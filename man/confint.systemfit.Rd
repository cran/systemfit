\name{confint.systemfit}
\alias{confint.systemfit}
\title{Confidence intervals of coefficients}

\description{
   This function calculates the confidence intervals of the
   coefficients of one equation
   from an object returned by \code{\link{systemfit}}.
}

\usage{
   \method{confint}{systemfit}( object, parm = NULL, level = 0.95, \dots )
}

\arguments{
  \item{object}{an object of type \code{systemfit}.}
  \item{parm}{not use yet.}
  \item{level}{confidence level.}
  \item{\dots}{other arguments.}
}

\value{
   An object of class \code{confint.systemfit}, which is a matrix
   with columns giving lower and upper confidence limits
   for each estimated coefficient. These will be labelled as
   (1-level)/2 and 1 - (1-level)/2 in % (by default 2.5% and 97.5%).
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{
   \code{\link{systemfit}}, \code{\link{confint.systemfit.equation}},
   \code{\link{print.confint.systemfit}} and \code{\link{confint}}
}


\examples{
\dontrun{library( systemfit )}

data( kmenta )
attach( kmenta )
demand <- q ~ p + d
supply <- q ~ p + f + a
labels <- list( "demand", "supply" )
system <- list( demand, supply )

## perform OLS on each of the equations in the system
fitols <- systemfit( "OLS", system, labels )

## print the confidence interval of the coefficients
confint( fitols )
}

\keyword{models}
\keyword{regression}



