
%  $Id: lrtest.systemfit.Rd 986 2008-10-16 07:57:24Z arne $

\name{lrtest.systemfit}
\alias{lrtest.systemfit}

\title{Likelihood Ratio test for Equation Systems}

\description{
   Testing linear hypothesis on the coefficients of a system of equations
   by a Likelihood Ratio test.
}

\usage{
   \method{lrtest}{systemfit}( object, ... )
}

\arguments{
   \item{object}{a fitted model object of class \code{systemfit}.}
   \item{\dots}{further fitted model objects of class \code{systemfit}.}
}

\details{
   \code{lrtest.systemfit} consecutively compares
   the fitted model object \code{object}
   with the models passed in \code{...}.

   The LR-statistic for sytems of equations is
   \deqn{
   LR = T \cdot \left(
      log \left| \hat{ \hat{ \Sigma } }_r \right|
      - log \left| \hat{ \hat{ \Sigma } }_u \right|
      \right)
   }
   where \eqn{T} is the number of observations per equation, and
   \eqn{\hat{\hat{\Sigma}}_r} and \eqn{\hat{\hat{\Sigma}}_u} are
   the residual covariance matrices calculated by formula "0"
   (see \code{\link{systemfit}})
   of the restricted and unrestricted estimation, respectively.
   Asymptotically, \eqn{LR} has a \eqn{\chi^2}
   distribution with \eqn{j} degrees of freedom
   under the null hypothesis
   (Green, 2003, p. 349).
}

\value{
   An object of class \code{anova},
   which contains the log-likelihood value,
   degrees of freedom, the difference in degrees of freedom,
   likelihood ratio Chi-squared statistic and corresponding p value.
   See documentation of \code{\link[lmtest]{lrtest}}
   in package "lmtest".
}

\references{
  Greene, W. H. (2003)
  \emph{Econometric Analysis, Fifth Edition}, Prentice Hall.
}

\author{Arne Henningsen \email{arne.henningsen@googlemail.com}}

\seealso{\code{\link{systemfit}}, \code{\link[lmtest]{lrtest}}
   (package "lmtest"),
   \code{\link{linear.hypothesis.systemfit}}}

\examples{
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
system <- list( demand = eqDemand, supply = eqSupply )

## unconstrained SUR estimation
fitsur <- systemfit( system, "SUR", data = Kmenta )

# create restriction matrix to impose \eqn{beta_2 = \beta_6}
R1 <- matrix( 0, nrow = 1, ncol = 7 )
R1[ 1, 2 ] <- 1
R1[ 1, 6 ] <- -1

## constrained SUR estimation
fitsur1 <- systemfit( system, "SUR", data = Kmenta, restrict.matrix = R1 )

## perform LR-test
lrTest1 <- lrtest( fitsur1, fitsur )
print( lrTest1 )   # rejected

# create restriction matrix to impose \eqn{beta_2 = - \beta_6}
R2 <- matrix( 0, nrow = 1, ncol = 7 )
R2[ 1, 2 ] <- 1
R2[ 1, 6 ] <- 1

## constrained SUR estimation
fitsur2 <- systemfit( system, "SUR", data = Kmenta, restrict.matrix = R2 )

## perform LR-test
lrTest2 <- lrtest( fitsur2, fitsur )
print( lrTest2 )   # accepted
}

\keyword{models}



