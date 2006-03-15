
%  $Id: lrtest.systemfit.Rd 171 2006-02-20 07:45:38Z henningsena $

\name{lrtest.systemfit}
\alias{lrtest.systemfit}
\alias{print.lrtest.systemfit}

\title{Likelihood Ratio test for Equation Systems}

\description{
   Likelihood Ratio test for linear parameter restrictions
   in equation system.
}

\usage{
   lrtest.systemfit( resultc, resultu )

   \method{print}{lrtest.systemfit}( x, digits = 4, ... )
}

\arguments{
   \item{resultc}{an object of type \code{systemfit} that
      contains the results of the restricted estimation.}
   \item{resultu}{an object of type \code{systemfit} that
      contains the results of the unconstrained estimation.}
   \item{x}{an object of class \code{ftest.systemfit}.}
   \item{digits}{number of digits to print.}
   \item{\dots}{currently not used.}
}

\details{
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
   \code{lrtest.systemfit} returns a list of class \code{lrtest.systemfit}
   that includes following objects:

   \item{statistic}{the empirical likelihood ratio statistic.}
   \item{p.value}{the p-value of the \eqn{\chi^2}-test.}
   \item{nRestr}{number of restrictions (\eqn{j}, degrees of freedom).}
}

\references{
  Greene, W. H. (2003)
  \emph{Econometric Analysis, Fifth Edition}, Prentice Hall.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{\code{\link{systemfit}}, \code{\link{ftest.systemfit}},
   \code{\link{waldtest.systemfit}}}

\examples{
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
system <- list( demand = eqDemand, supply = eqSupply )

## unconstrained SUR estimation
fitsur <- systemfit( "SUR", system, data = Kmenta )

# create restriction matrix to impose \eqn{beta_2 = \beta_6}
R1 <- matrix( 0, nrow = 1, ncol = 7 )
R1[ 1, 2 ] <- 1
R1[ 1, 6 ] <- -1

## constrained SUR estimation
fitsur1 <- systemfit( "SUR", system, data = Kmenta, R.restr = R1 )

## perform LR-test
lrTest1 <- lrtest.systemfit( fitsur1, fitsur )
print( lrTest1 )   # rejected

# create restriction matrix to impose \eqn{beta_2 = - \beta_6}
R2 <- matrix( 0, nrow = 1, ncol = 7 )
R2[ 1, 2 ] <- 1
R2[ 1, 6 ] <- 1

## constrained SUR estimation
fitsur2 <- systemfit( "SUR", system, data = Kmenta, R.restr = R2 )

## perform LR-test
lrTest2 <- lrtest.systemfit( fitsur2, fitsur )
print( lrTest2 )   # accepted
}

\keyword{models}



