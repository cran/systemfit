\name{print.hausman.systemfit}
\alias{print.hausman.systemfit}
\title{Print result of Hausman test}

\description{
   This function prints the results of a Hausman test.
}

\usage{
   \method{print}{hausman.systemfit}( x, digits=6, \dots )
}

\arguments{
   \item{x}{an object of type \code{hausman.systemfit}.}
   \item{digits}{number of digits to print.}
   \item{\dots}{other arguments.}
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de} }

\seealso{\code{\link{hausman.systemfit}}}


\examples{
\dontrun{library( systemfit )}

data( kmenta )
demand <- q ~ p + d
supply <- q ~ p + f + a
inst <- ~ d + f + a
labels <- list( "demand", "supply" )
system <- list( demand, supply )

## perform the estimations
fit2sls <- systemfit( "2SLS", system, labels, inst, data = kmenta )
fit3sls <- systemfit( "3SLS", system, labels, inst, data = kmenta )

## perform the Hausman test
h <- hausman.systemfit( fit2sls, fit3sls )
print( h )
}

\keyword{models}



