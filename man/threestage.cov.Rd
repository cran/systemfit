
% 	$Id: threestage.cov.Rd,v 1.4 2002/12/18 06:40:09 hamannj Exp $	

\name{threestage.cov}
\alias{threestage.cov}
\alias{threestage.cov.systemfit}

\title{Variance-Covariance Matrix}

\description{
  This function returns a variance-covariance estimation matrix from a
  resulting simultaneous estimation object such as type
  \code{threestage}.
}

\usage{
threestage.cov( results, eqni, eqnj )
}

\arguments{
  \item{results}{a set of 3SLS objects returned from \code{threestage} }
  \item{eqni}{the ith equation in the set of \code{threestage} objects}
  \item{eqnj}{the jth equation in the set of \code{threestage} objects}
}

\value{
  \code{threestage.cov} returns a submatrix from the variance-covariance
  matrix from the variance-covariance matrix used for estimation during
  3SLS.}
}

\references{

  Hasenauer, H; Monserud, R and T. Gregoire. (1998)
  Using Simultansous Regression Techniques with Individual-Tree Growth Models.
  \emph{Forest Science}. 44(1):87-95
  
}

\author{Jeff D. Hamann \email{jeff\_hamann@hamanndonald.com}}

\seealso{\code{\link{ols}},\code{\link{twostage}},\code{\link{sur}}and \code{\link{threestage}}}


\examples{
library( systemfit )

data( kmenta )
attach( kmenta )
demand <- q ~ p + d
supply <- q ~ p + f + a
inst <- ~ d + f + a
labels <- list( "demand", "supply" )
system <- list( demand, supply )

## perform the estimation and report the results for the whoel system
fit3sls <- threestage.systemfit( system, inst, labels, kmenta )
print( fit3sls )

## get the variance-covariance matrix used for estimation 
print( "covariance of residuals used for estimation (from 2sls)" )
print( threestage.cov( fit3sls, 1, 2 ) )

}
\keyword{models}



