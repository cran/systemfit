
% 	$Id: lrtest.systemfit.Rd,v 1.3 2004/02/23 13:06:55 henningsena Exp $	

\name{lrtest.systemfit}
\alias{lrtest.systemfit}

\title{Likelihood Ratio test for Equation Systems}

\description{
  Perform a Likelihood Ratio test to test linear restrictions
  when estimating an equation system.
}

\usage{
lrtest.systemfit( resultc, resultu )
}

\arguments{
  \item{resultc}{an object of type \code{systemfit} that
    contains the results of the restricted estimation.}
  \item{resultu}{an object of type \code{systemfit} that
    contains the results of the unconstrained estimation.}
}

\details{
  At the moment only the LR test for SUR models is implemented.

  This test calculates the likelihood ratio value by calculating
  the estimated residual covariance matrix \eqn{\hat{\Sigma}} by
  the maximum likelihood formula
  \eqn{\hat{\sigma}_{ij} = (\hat{e}_i' \hat{e}_j) / T}
  (no matter which formula is used in the estimation).
}

\value{
  \code{lrtest.sur.systemfit} returns a three values:

  \item{lr}{the empirical likelihood ratio value.}
  \item{p}{the p-value of a \eqn{\chi^2}-test.}
  \item{df}{the degrees of freedom of the test (= number of restrictions).}
}

\references{

  Greene, W. H. (2002)
  \emph{Econometric Analysis, Fifth Edition}, Prentice Hall.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{\code{\link{systemfit}}}

\examples{
\dontrun{library( systemfit )}
data( kmenta )
demand <- q ~ p + d
supply <- q ~ p + f + a
labels <- list( "demand", "supply" )
system <- list( demand, supply )

## unconstrained SUR estimation
fitsur <- systemfit("SUR", system, labels, data=kmenta )
print( fitsur )

## SUR estimation with 2 restrictions
Rrestr <- matrix(0,2,7)
qrestr <- matrix(0,2,1)
Rrestr[1,3] <-  1
Rrestr[1,7] <- -1
Rrestr[2,2] <- -1
Rrestr[2,5] <-  1
qrestr[2,1] <-  0.5
fitsur2 <- systemfit("SUR", system, labels, data=kmenta,
                      R.restr=Rrestr, q.restr=qrestr )
print( fitsur2 )

## perform LR-test
lrtest.systemfit( fitsur2, fitsur )
}

\keyword{models}



