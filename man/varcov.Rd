
% 	$Id: varcov.Rd,v 1.3 2002/12/18 06:40:09 hamannj Exp $	

\name{varcov}
\alias{varcov}
\alias{varcov.systemfit}

\title{Variance-Covariance}

\description{
  The function returns the variance-covaraince of the residuals for a
  set of equations from the residuals. The values of the elements are
  defined as,

  \deqn{\frac{e_{i}'e_{j}}{\sqrt{n_i n_j}}}
  
  where $e_i$ and $e_j$ are the residuals and $n_i$ and $n_j$ are the
  error degrees of freedom for equations i and j.
}

\usage{
varcov.systemfit( results )
}

\arguments{
  \item{results}{an object of type \code{twostage.systemfit}.}
}

\value{
  \code{varcov} returns a variance-covaraince matrix of the residuals
  from a set of objects of class \code{twostage}.
}

\references{    

  Greene, W. H. (1993)
  \emph{Econometric Analysis, Second Edition}, Macmillan.
  
  Kmenta, J. (1997)
  \emph{Elements of Econometrics, Second Edition}, University of
  Michigan Publishing

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

## perform 2SLS on each of the equations in the system
fit2sls <- twostage.systemfit( system, inst, labels, kmenta )
fit3sls <- threestage.systemfit( system, inst, labels, kmenta )

## print the results
print( fit2sls )
print( fit3sls )
print( "covariance of residuals used for estimation (from 2sls)" )
print( varcov.systemfit( fit2sls ) )
print( "covariance of residuals" )
print( varcov.systemfit( fit3sls ) )

}
\keyword{models}



