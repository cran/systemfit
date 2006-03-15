\name{ftest.systemfit}
\alias{ftest.systemfit}
\alias{print.ftest.systemfit}

\title{F-test for Equation Systems}

\description{
   F-test for linear parameter restrictions
   in equation systems.
}

\usage{
   ftest.systemfit( object, R.restr,
      q.restr = rep( 0, nrow( R.restr ) ) )

   \method{print}{ftest.systemfit}( x, digits = 4, ... )
}

\arguments{
   \item{object}{an object of type \code{systemfit}.}
   \item{R.restr}{j x k matrix to impose linear
      restrictions on the parameters by \code{R.restr} * \eqn{b} = \code{q.restr}
      (j = number of restrictions, k = number of all parameters,
      \eqn{b} = vector of all parameters).}
   \item{q.restr}{an optional vector with j elements to impose linear
      restrictions (see \code{R.restr}); default is a vector
      that contains only zeros.}
   \item{x}{an object of class \code{ftest.systemfit}.}
   \item{digits}{number of digits to print.}
   \item{\dots}{currently not used.}
}

\details{
   The F-statistic for sytems of equations is
   \deqn{F = \frac{
      ( R \hat{b} - q )'
      ( R ( X' ( \hat{\Sigma} \otimes I )^{-1} X )^{-1} R' )^{-1}
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
   \eqn{\Sigma} is the residual covariance matrix used in the estimation.
   Under the null hypothesis, \eqn{F} has an F-distribution
   with \eqn{j} and \eqn{M \cdot T - K} degrees of freedom
   (Theil, 1971, p. 314).
}

\value{
   \code{ftest.systemfit} returns a list of class \code{ftest.systemfit}
   that includes following objects:

   \item{statistic}{the empirical F statistic.}
   \item{p.value}{the p-value of the F-test.}
   \item{nRestr}{number of restrictions (\eqn{j}, degrees of freedom of the
      numerator).}
   \item{dfSys}{degrees of freedom of the equation system
      (\eqn{M \cdot T - K}, degrees of freedom of the denominator).}
}

\references{
Theil, Henri (1971). \emph{Principles of Econometrics}, John Wiley & Sons, New York.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{\code{\link{systemfit}}, \code{\link{waldtest.systemfit}},
   \code{\link{lrtest.systemfit}}}

\examples{
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
system <- list( demand = eqDemand, supply = eqSupply )

## unconstrained SUR estimation
fitsur <- systemfit( "SUR", system, data=Kmenta )

# create restriction matrix to test whether \eqn{beta_2 = \beta_6}
R1 <- matrix( 0, nrow = 1, ncol = 7 )
R1[ 1, 2 ] <- 1
R1[ 1, 6 ] <- -1

## perform F-test
fTest1 <- ftest.systemfit( fitsur, R1 )
print( fTest1 )   # rejected

# create restriction matrix to test whether \eqn{beta_2 = - \beta_6}
R2 <- matrix( 0, nrow = 1, ncol = 7 )
R2[ 1, 2 ] <- 1
R2[ 1, 6 ] <- 1

## perform F-test
fTest2 <- ftest.systemfit( fitsur, R2 )
print( fTest2 )   # accepted
}

\keyword{models}



