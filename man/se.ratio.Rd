\name{se.ratio}
\alias{se.ratio}
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
  \item{resultsi}{an object of type systemfit.system (\code{ols}, \code{twostage} or \code{threestage}.}
  \item{resultsj}{an object of type systemfit.system (\code{ols}, \code{twostage} or \code{threestage}.}
  \item{eqni}{index for equation to obtain the ratio of standard errors}
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
fit2sls <- twostage.systemfit( system, inst, labels, kmenta )
fit3sls <- threestage.systemfit( system, inst, labels, kmenta )

## print the results from the fits
print( fit2sls )
print( fit3sls )
print( "covariance of residuals used for estimation (from 2sls)" )
print( varcov.systemfit( fit2sls ) )
print( "covariance of residuals" )
print( varcov.systemfit( fit3sls ) )

## examine the correlation between the predicted values
## of suppy and demand by plotting the correlation over
## the value of q
r12 <- correlation.systemfit( fit3sls, 1, 2 )
plot( q, r12, main="correlation between predictions from supply and demand" )  

## examine the improvement of 3SLS over OLS by computing
## the ratio of the standard errors of the estimates
improve.ratio <- se.ratio.systemfit( fit2sls, fit3sls, 2 )
print( "summary values for the ratio in the std. err. for the predictions" )
print( summary( improve.ratio ) )
}
\keyword{models}



