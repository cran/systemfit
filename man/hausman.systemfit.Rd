
%  $Id: hausman.systemfit.Rd 986 2008-10-16 07:57:24Z arne $

\name{hausman.systemfit}
\alias{hausman.systemfit}

\title{Hausman Test}

\description{
  \code{hausman.systemfit} returns the Hausman statistic for a specification test.
}

\usage{
   hausman.systemfit( results2sls, results3sls )
}

\arguments{
  \item{ results2sls }{result of a \emph{2SLS} (limited information) estimation
     returned by \code{\link{systemfit}}.}
  \item{ results3sls }{result of a \emph{3SLS} (full information) estimation
     returned by \code{\link{systemfit}}.}
}

\details{
   The null hypotheses of the test is that all exogenous variables are
   uncorrelated with all disturbance terms.
   Under this hypothesis both the 2SLS and the 3SLS estimator are consistent
   but only the 3SLS estimator is (asymptotically) efficient.
   Under the alternative hypothesis the 2SLS estimator is consistent
   but the 3SLS estimator is inconsistent.

   The Hausman test statistic is

   \deqn{m = ( b_2 - b_3 )' ( V_2 - V_3 ) ( b_2 - b_3 ) }

   where $b_2$ and $V_2$ are the estimated coefficients and their variance
   covariance matrix of a \emph{2SLS} estimation and
   $b_3$ and $V_3$ are the estimated coefficients and their variance
   covariance matrix of a \emph{3SLS} estimation.
}

\value{
   \code{hausman.systemfit} returns a list of the class
    \code{htest} that contains following elements:
   \item{q}{vector of the differences between the estimated coefficients.}
   \item{qVar}{variance covariance matrix of \code{q} (difference between the
     variance covariance matrices of the estimated coefficients).}
   \item{statistic}{the Hausman test statistic.}
   \item{parameter}{degrees of freedom.}
   \item{p.value}{P-value of the test.}
   \item{method}{character string describing this test.}
   \item{data.name}{name of the data.frame used for estimation.}
}

\references{

  Greene, W. H. (1993)
  \emph{Econometric Analysis, Fifth Edition}, Macmillan.

  Hausman, J. A. (1978)
  Specification Tests in Econometrics. \emph{Econometrica}. 46:1251-1271.

  Kmenta, J. (1997)
  \emph{Elements of Econometrics, Second Edition}, University of
  Michigan Publishing

}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com},\cr
  Arne Henningsen \email{arne.henningsen@googlemail.com}
}

\seealso{\code{\link{systemfit}}}

\examples{
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
inst <- ~ income + farmPrice + trend
system <- list( demand = eqDemand, supply = eqSupply )

## perform the estimations
fit2sls <- systemfit( system, "2SLS", inst = inst, data = Kmenta )
fit3sls <- systemfit( system, "3SLS", inst = inst, data = Kmenta )

## perform the Hausman test
h <- hausman.systemfit( fit2sls, fit3sls )
print( h )
}

\keyword{models}



