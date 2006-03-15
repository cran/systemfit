\name{waldtest.systemfit}
\alias{waldtest.systemfit}
\alias{print.waldtest.systemfit}

\title{Wald-test for Equation Systems}

\description{
   Wald-test for linear parameter restrictions
   in equation systems.
}

\usage{
   waldtest.systemfit( object, R.restr,
      q.restr = rep( 0, nrow( R.restr ) ) )

   \method{print}{waldtest.systemfit}( x, digits = 4, ... )
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
   \item{x}{an object of class \code{waldtest.systemfit}.}
   \item{digits}{number of digits to print.}
   \item{\dots}{currently not used.}
}

\details{
   The Wald-statistic for sytems of equations is
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
   \code{waldtest.systemfit} returns a list of class \code{waldtest.systemfit}
   that includes following objects:

   \item{statistic}{the empirical Wald statistic.}
   \item{p.value}{the p-value of the Wald-test.}
   \item{nRestr}{number of restrictions (\eqn{j}, degrees of freedom).}
}

\references{
  Greene, W. H. (2003)
  \emph{Econometric Analysis, Fifth Edition}, Prentice Hall.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{\code{\link{systemfit}}, \code{\link{ftest.systemfit}},
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

## perform Wald-test
waldTest1 <- waldtest.systemfit( fitsur, R1 )
print( waldTest1 )   # rejected

# create restriction matrix to test whether \eqn{beta_2 = - \beta_6}
R2 <- matrix( 0, nrow = 1, ncol = 7 )
R2[ 1, 2 ] <- 1
R2[ 1, 6 ] <- 1

## perform Wald-test
waldTest2 <- waldtest.systemfit( fitsur, R2 )
print( waldTest2 )   # accepted
}

\keyword{models}



