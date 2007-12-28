
%  $Id: se.ratio.systemfit.Rd 437 2007-06-30 07:44:12Z henningsena $


\name{se.ratio.systemfit}
\alias{se.ratio.systemfit}

\title{Ratio of the Standard Errors}

\description{
  \code{se.ratio.systemfit} returns a vector of the ratios of the
  standard errors of the predictions for two equations.
}

\usage{
se.ratio.systemfit( resultsi, resultsj, eqni )
}

\arguments{
  \item{resultsi}{an object of type \code{systemfit}.}
  \item{resultsj}{an object of type \code{systemfit}.}
  \item{eqni}{index for equation to obtain the ratio of standard errors.}
}

\value{
  \code{se.ratio} returns a vector of the standard errors of the ratios
  for the predictions between the predicted values in equation i and
  equation j.

}

\references{

  Hasenauer, H; Monserud, R and T. Gregoire. (1998)
  Using Simultaneous Regression Techniques with Individual-Tree Growth Models.
  \emph{Forest Science}. 44(1):87-95

}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com}}

\seealso{\code{\link{systemfit}} and \code{\link{correlation.systemfit}}}


\examples{
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
inst <- ~ income + farmPrice + trend
system <- list( demand = eqDemand, supply = eqSupply )

## perform 2SLS on each of the equations in the system
fit2sls <- systemfit( system, "2SLS", inst = inst, data = Kmenta )
fit3sls <- systemfit( system, "3SLS", inst = inst, data = Kmenta )

## print the results from the fits
print( fit2sls )
print( fit3sls )
print( "covariance of residuals used for estimation (from 2sls)" )
print( fit3sls$residCovEst )
print( "covariance of residuals" )
print( fit3sls$residCov )

## examine the improvement of 3SLS over 2SLS by computing
## the ratio of the standard errors of the estimates
improve.ratio <- se.ratio.systemfit( fit2sls, fit3sls, 2 )
print( "summary values for the ratio in the std. err. for the predictions" )
print( summary( improve.ratio ) )
}
\keyword{models}



