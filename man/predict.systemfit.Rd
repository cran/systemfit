\name{predict.systemfit}
\alias{predict.systemfit}
\alias{predict.systemfit.equation}

\title{Predictions from System Estimation}

\description{
   Returns the predicted values, their standard errors and the
   confidence limits of prediction.
}

\usage{
\method{predict}{systemfit}( object, newdata = NULL,
                             se.fit = FALSE, se.pred = FALSE,
                             interval = "none", level=0.95,
                             useDfSys = NULL, ... )

\method{predict}{systemfit.equation}( object, newdata = NULL,
                             se.fit = FALSE, se.pred = FALSE,
                             interval = "none", level=0.95,
                             useDfSys = NULL, ... )
}

\arguments{
   \item{object}{an object of class \code{systemfit} or
      \code{systemfit.equation}.}
   \item{newdata}{An optional data frame in which to look for variables with
      which to predict.  If it is \code{NULL}, the fitted values are returned.}
   \item{se.fit}{return the standard error of the fitted values?}
   \item{se.pred}{return the standard error of prediction?}
   \item{interval}{Type of interval calculation ("none", "confidence"
      or "prediction")}
   \item{level}{Tolerance/confidence level.}
   \item{useDfSys}{logical. Use the degrees of freedom of the whole system
      (in place of the degrees of freedom of the single equation)
      to calculate the confidence or prediction intervals.
      If it not specified (\code{NULL}), it is set to \code{TRUE}
      if restrictions on the coefficients are imposed and \code{FALSE} otherwise.}
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
   ("<eqnLable>.pred") and if requested
   the standard errors of the fitted values ("<eqnLable>.se.fit"),
   the standard errors of the prediction ("<eqnLable>.se.pred"),
   and the lower ("<eqnLable>.lwr") and upper ("<eqnLable>.upr")
   limits of the confidence or prediction interval(s).

   \code{predict.systemfit.equation} returns a dataframe that
   contains the predicted values ("fit") and if requested
   the standard errors of the fitted values ("se.fit"),
   the standard errors of the prediction ("se.pred"),
   and the lower ("lwr") and upper ("upr")
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

\author{Arne Henningsen \email{arne.henningsen@googlemail.com}}


\seealso{
   \code{\link{systemfit}}, \code{\link{predict}}
}

\examples{
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
system <- list( demand = eqDemand, supply = eqSupply )

## OLS estimation
fitols <- systemfit( system, data=Kmenta )

## predicted values and limits
predict( fitols )

## predicted values of the first equation
predict( fitols$eq[[1]] )

## predicted values of the second equation
predict( fitols$eq[[2]] )
}

\keyword{models}



