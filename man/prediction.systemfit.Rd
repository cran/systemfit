
% 	$Id: prediction.systemfit.Rd,v 1.3 2003/11/18 08:54:16 henningsena Exp $	

\name{prediction.systemfit}
\alias{prediction.systemfit}

\title{Predictions from Equation System Estimation}

\description{
   Returns the predicted values, their standard errors and the
   confidence limits of prediction.
}

\usage{
prediction.systemfit( object, data=object$data, alpha=0.05 )
}

\arguments{
  \item{object}{an object of type systemfit.system.}
  \item{data}{a dataframe containing the exogenous variables.}
  \item{alpha}{level of significance.}
}

\value{
  \code{prediction.systemfit} returns a list of objects for each
  equation that contains the following components:

  \item{predicted}{predicted values.}
  \item{se.prediction}{standard error of the predicted values.}
  \item{prediction.limits}{lower and upper (1-alpha) confidence limits of prediction.}
}

\references{

  Greene, W. H. (1993)
  \emph{Econometric Analysis, Second Edition}, Macmillan.

  Kmenta, J. (1997)
  \emph{Elements of Econometrics, Second Edition}, University of
  Michigan Publishing.
}

\author{Arne Henningsen (\url{http://www.arne-henningsen.de})}


}

\seealso{\code{\link{systemfit}}}

\examples{
\dontrun{library( systemfit )}
data( kmenta )
demand <- q ~ p + d
supply <- q ~ p + f + a
labels <- list( "demand", "supply" )
system <- list( demand, supply )

## OLS estimation
fitols <- systemfit("OLS", system, labels, data=kmenta )
print( fitols )

## calculate predicted values and limits
pred <- prediction.systemfit( fitols )
print( pred[[2]]$prediction.limits )
}

\keyword{models}



