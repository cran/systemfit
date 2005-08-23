
% 	$Id: correlation.systemfit.Rd,v 1.4 2004/02/23 13:06:55 henningsena Exp $	

\name{correlation.systemfit}
\alias{correlation.systemfit}

\title{Correlation between Predictions from Equation i and j}

\description{
  \code{correlation} returns a vector of the correlations
  between the preditions of two equations in a set of equations. The
  correlation between the predictions is defined as,

  \deqn{
    r_{ijk} = \frac{x'_{ik}C_{ij}x_{jk}}{\sqrt{(x'_{ik}C_{ii}x_{ik})(x'_{jk}C_{jj}x_{jk})}}
  }

  where \eqn{r_{ijk}} is the correlation between the predicted values of
  equation i and j and  \eqn{C_{ij}} is the cross-equation variance-covariance
  matrix between equations i and j.

}

\usage{
correlation.systemfit( results, eqni, eqnj )
}

\arguments{
  \item{results}{an object of type \code{systemfit}.}
  \item{eqni}{index for equation i}
  \item{eqnj}{index for equation j}
}

\value{
  \code{correlation} returns a vector of the correlations between the
  predicted values in equation i and equation j.
}

\references{    

  Greene, W. H. (1993)
  \emph{Econometric Analysis, Second Edition}, Macmillan.

  Hasenauer, H; Monserud, R and T. Gregoire. (1998)
  Using Simultansous Regression Techniques with Individual-Tree Growth Models.
  \emph{Forest Science}. 44(1):87-95

  Kmenta, J. (1997)
  \emph{Elements of Econometrics, Second Edition}, University of
  Michigan Publishing

}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com}}

\seealso{\code{\link{systemfit}}}


\examples{
\dontrun{library( systemfit )}

data( kmenta )
attach( kmenta )
demand <- q ~ p + d
supply <- q ~ p + f + a
inst <- ~ d + f + a
labels <- list( "demand", "supply" )
system <- list( demand, supply )

## perform 2SLS on each of the equations in the system
fit2sls <- systemfit( "2SLS", system, labels, inst)
print( fit2sls )
print( fit2sls$rcov )

## perform the 3SLS
fit3sls <- systemfit( "3SLS", system, labels, inst )
print( fit3sls )
print( "covariance of residuals used for estimation (from 2sls)" )
print( fit3sls$rcovest )
print( "covariance of residuals" )
print( fit3sls$rcov )

## examine the correlation between the predicted values
## of suppy and demand by plotting the correlation over
## the value of q
r12 <- correlation.systemfit( fit3sls, 1, 2 )
plot( q, r12, main="correlation between predictions from supply and demand" )
}

\keyword{models}
