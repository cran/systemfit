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

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com}}

\seealso{\code{\link{systemfit}}, \code{\link{print.systemfit.system}} and
  \code{\link{print.systemfit.system}}}


\examples{
\dontrun{library( systemfit )}

data( kmenta )
attach( kmenta )
demand <- q ~ p + d
supply <- q ~ p + f + a
inst <- ~ d + f + a
labels <- list( "demand", "supply" )
system <- list( demand, supply )

## perform OLS on each of the equations in the system
fitols <- systemfit( "OLS", system, labels )

## print the results
summary.systemfit.system( fitols )
}

\keyword{models}



