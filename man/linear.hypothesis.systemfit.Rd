\name{linear.hypothesis.systemfit}
\alias{linear.hypothesis.systemfit}

\title{Test Linear Hypothesis}

\description{
   Testing linear hypothesis on parameters of a system of equations
   by an F-test or Wald-test.
}

\usage{
   \method{linear.hypothesis}{systemfit}( model,
      hypothesis.matrix, rhs = NULL, test = c( "FT", "F", "Chisq" ),
      vcov. = NULL, ... )
}

\arguments{
   \item{model}{a fitted object of type \code{systemfit}.}
   \item{hypothesis.matrix}{matrix (or vector) giving linear combinations
      of coefficients by rows,
      or a character vector giving the hypothesis in symbolic form
      (see documentation of \code{\link[car]{linear.hypothesis}}
      in package "car" for details).}
   \item{rhs}{optional right-hand-side vector for hypothesis,
      with as many entries as rows in the hypothesis matrix;
      if omitted, it defaults to a vector of zeroes.}
   \item{test}{character string, "\code{FT}", "\code{F}", or "\code{Chisq}",
      specifying whether to compute
      Theil's finite-sample F test (with approximate F distribution),
      the finite-sample Wald test (with approximate F distribution),
      or the large-sample Wald test
      (with asymptotic Chi-squared distribution).}
   \item{vcov.}{a function for estimating the covariance matrix
      of the regression coefficients or an estimated covariance matrix
      (function \code{vcov} is used by default).}
   \item{\dots}{further arguments passed to
      \code{\link[car]{linear.hypothesis.default}} (package "car").}
}

\details{
   Theil's \eqn{F} statistic for sytems of equations is
   \deqn{F = \frac{
      ( R \hat{b} - q )'
      ( R ( X' ( \Sigma \otimes I )^{-1} X )^{-1} R' )^{-1}
      ( R \hat{b} - q ) /
      j
   }{
      \hat{e}' ( \Sigma \otimes I )^{-1} \hat{e} /
      ( M \cdot T - K )
   } }
   where \eqn{j} is the number of restrictions,
   \eqn{M} is the number of equations,
   \eqn{T} is the number of observations per equation,
   \eqn{K} is the total number of estimated coefficients, and
   \eqn{\Sigma} is the estimated residual covariance matrix.
   Under the null hypothesis, \eqn{F} has an approximate \eqn{F} distribution
   with \eqn{j} and \eqn{M \cdot T - K} degrees of freedom
   (Theil, 1971, p. 314).

   The \eqn{F} statistic for a Wald test is
   \deqn{ F = \frac{
      ( R \hat{b} - q )'
      ( R \, \widehat{Cov} [ \hat{b} ] R' )^{-1}
      ( R \hat{b} - q )
   }{
      j
   } }
   Under the null hypothesis, \eqn{F} has an approximate \eqn{F} distribution
   with \eqn{j} and \eqn{M \cdot T - K} degrees of freedom
   (Greene, 2003, p. 346).

   The \eqn{\chi^2} statistic for a Wald test is
   \deqn{ W =
      ( R \hat{b} - q )'
      ( R \widehat{Cov} [ \hat{b} ] R' )^{-1}
      ( R \hat{b} - q )
   }
   Asymptotically, \eqn{W} has a \eqn{\chi^2}
   distribution with \eqn{j} degrees of freedom
   under the null hypothesis
   (Greene, 2003, p. 347).
}

\value{
   An object of class \code{anova},
   which contains the residual degrees of freedom in the model,
   the difference in degrees of freedom,
   the test statistic (either F or Wald/Chisq)
   and the corresponding p value.
   See documentation of \code{\link[car]{linear.hypothesis}}
   in package "car".
}

\references{
  Greene, W. H. (2003)
  \emph{Econometric Analysis, Fifth Edition}, Prentice Hall.

   Theil, Henri (1971)
   \emph{Principles of Econometrics}, John Wiley & Sons, New York.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{\code{\link{systemfit}}, \code{\link[car]{linear.hypothesis}}
   (package "car"),
   \code{\link{lrtest.systemfit}}}

\examples{
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
system <- list( demand = eqDemand, supply = eqSupply )

## unconstrained SUR estimation
fitsur <- systemfit( system, method = "SUR", data=Kmenta )

# create hypothesis matrix to test whether beta_2 = \beta_6
R1 <- matrix( 0, nrow = 1, ncol = 7 )
R1[ 1, 2 ] <- 1
R1[ 1, 6 ] <- -1
# the same hypothesis in symbolic form
restrict1 <- "demand_price - supply_farmPrice = 0"

## perform Theil's F test
linear.hypothesis( fitsur, R1 )  # rejected
linear.hypothesis( fitsur, restrict1 )

## perform Wald test with F statistic
linear.hypothesis( fitsur, R1, test = "F" )  # rejected
linear.hypothesis( fitsur, restrict1 )

## perform Wald-test with chi^2 statistic
linear.hypothesis( fitsur, R1, test = "Chisq" )  # rejected
linear.hypothesis( fitsur, restrict1, test = "Chisq" )

# create hypothesis matrix to test whether beta_2 = - \beta_6
R2 <- matrix( 0, nrow = 1, ncol = 7 )
R2[ 1, 2 ] <- 1
R2[ 1, 6 ] <- 1
# the same hypothesis in symbolic form
restrict2 <- "demand_price + supply_farmPrice = 0"

## perform Theil's F test
linear.hypothesis( fitsur, R2 )  # accepted
linear.hypothesis( fitsur, restrict2 )

## perform Wald test with F statistic
linear.hypothesis( fitsur, R2, test = "F" )  # accepted
linear.hypothesis( fitsur, restrict2 )

## perform Wald-test with chi^2 statistic
linear.hypothesis( fitsur, R2, test = "Chisq" )  # accepted
linear.hypothesis( fitsur, restrict2, test = "Chisq" )
}

\keyword{models}



