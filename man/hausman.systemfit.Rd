
% 	$Id: hausman.systemfit.Rd,v 1.5 2004/02/18 17:35:41 henningsena Exp $	

\name{hausman.systemfit}
\alias{hausman.systemfit}

\title{Hausman's Test}

\description{
  \code{hausman} returns the Hausman's statistic for specification.

  \deqn{q'(V_1 - V_0)q}

  where $V_1$ and $V_0$ are the covb values from a \code{twostage} or
  \code{threestage} object and q is the difference in the b vectors from
  the \code{twostage} or \code{threestage} objects.

}

\usage{
hausman.systemfit( li.results, fi.results )
}

\arguments{
  \item{li.results}{the limited information results object
    (\code{twostage}) }
  \item{fi.results}{the full information results object
    (\code{threestage}) }
}

\value{
  \code{hausman.systemfit} returns the value of the test statistic.
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

## perform the estimation and report the results for the whole system
fit2sls <- systemfit( "2SLS", system, labels, inst)
fit3sls <- systemfit( "3SLS", system, labels, inst)

## perform the hausman test on the first equation
h <- hausman.systemfit( fit2sls, fit3sls )
pval <- 1 - pchisq( h, dim( fit3sls$bcov )[1] )

print( h )
print( pval )

}
\keyword{models}



