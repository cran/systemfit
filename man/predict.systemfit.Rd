\name{predict.systemfit}
\alias{predict.systemfit}

\title{Predictions from Equation System Estimation}

\description{
   Returns the predicted values, their standard errors and the
   confidence limits of prediction.
}

\usage{
\method{predict}{systemfit}( object, data = object$data,
                             se.fit = FALSE, se.pred = FALSE,
                             interval = "none", level=0.95, ... )
}

\arguments{
  \item{object}{an object of type \code{systemfit}.}
  \item{data}{data frame in which to predict.}
  \item{se.fit}{return the standard error of the fitted values?}
  \item{se.pred}{return the standard error of prediction?}
  \item{interval}{Type of interval calculation ("none", "confidence"
      or "prediction")}
  \item{level}{Tolerance/confidence level.}
  \item{...}{additional optional arguments.}
}

\details{
   The variance of the fitted values
   (used to calculate the standard errors of the fitted values
   and the "confidence interval") is calculated by
   \eqn{Var[E[y^0]-\hat{y}^0]=x^0 \; Var[b] \; {x^0}'}\cr
   The variances of the predicted values
   (used to calculate the standard errors of the predicted values
   and the "prediction intervals") is calculated by
   \eqn{Var[y^0-\hat{y}^0]=\hat{\sigma}^2+x^0 \; Var[b] \; {x^0}'}
}

\value{
   \code{predict.systemfit} returns a dataframe that
   contains for each equation the predicted values
   (e.g. "eq1.pred") and if requested
   the standard errors of the fitted values (e.g. "eq1.se.fit"),
   the standard errors of the prediction (e.g. "eq1.se.pred"),
   and the lower (e.g. "eq1.lwr") and upper (e.g. "eq1.upr")
   limits of the confidence or prediction interval(s).
}

\references{

  Greene, W. H. (2003)
  \emph{Econometric Analysis, Fifth Edition}, Macmillan.

  Gujarati, D. N. (1995)
  \emph{Basic Econometrics, Third Edition}, McGraw-Hill.

  Kmenta, J. (1997)
  \emph{Elements of Econometrics, Second Edition}, University of
  Michigan Publishing.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}


\seealso{
   \code{\link{systemfit}}, \code{\link{predict.systemfit.equation}}
   and \code{\link{predict}}
}

\examples{
\dontrun{library( systemfit )}
data( kmenta )
demand <- q ~ p + d
supply <- q ~ p + f + a
labels <- list( "demand", "supply" )
system <- list( demand, supply )

## OLS estimation
fitols <- systemfit("OLS", system, labels, data=kmenta )

## calculate predicted values and limits
predict( fitols )
}

\keyword{models}



