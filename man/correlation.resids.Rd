
% 	$Id: correlation.resids.Rd,v 1.1 2002/12/18 06:28:48 hamannj Exp $	


\name{correlation.resids}
\alias{correlation.resids}
\alias{correlation.resids.systemfit}

\title{Correlation matrix of the residuals}

\description{
  \code{correlation.resids.systemfit} returns a matrix of the
  correlations of the residuals for the systemfit object.
}

\usage{
correlation.resids.systemfit( results )
}

\arguments{
  \item{results}{an object of type systemfit.system (\code{ols},
    \code{twostage}, \code{sur}, or \code{threestage}.}
}

\value{
  \code{correlation.resids.systemfit} returns a matrix of the
  correlations of the residuals for the systemfit object.
}


%
%\references{    
%}
%

\author{Jeff D. Hamann \email{jeff\_hamann@hamanndonald.com}}

\seealso{\code{\link{ols}},\code{\link{twostage}}, \code{\link{sur}} and \code{\link{threestage}}}


\examples{
library( systemfit )

data( kmenta )
attach( kmenta )
demand <- q ~ p + d
supply <- q ~ p + f + a
inst <- ~ d + f + a
labels <- list( "demand", "supply" )
system <- list( demand, supply )

## perform 3SLS on each of the equations in the system
fit3sls <- threestage.systemfit( system, inst, labels, kmenta )

## print the results
print( fit3sls )

}
\keyword{models}



