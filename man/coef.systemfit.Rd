\name{coef.systemfit}
\alias{coef.systemfit}
\alias{coef.systemfit.equation}
\alias{coef.summary.systemfit}
\alias{coef.summary.systemfit.equation}
\title{Coefficients of systemfit object}

\description{
   These functions extract the coefficients from an object returned by
   \code{\link{systemfit}}.
}

\usage{
   \method{coef}{systemfit}( object, modified.regMat = FALSE, \dots )

   \method{coef}{systemfit.equation}( object, \dots )

   \method{coef}{summary.systemfit}( object, modified.regMat = FALSE, \dots )

   \method{coef}{summary.systemfit.equation}( object, \dots )
}

\arguments{
   \item{object}{an object of class \code{systemfit},
      \code{systemfit.equation}, \code{summary.systemfit}, or
      \code{summary.systemfit.equation}.}
   \item{modified.regMat}{logical. If \code{TRUE},
      the coefficients of the modified regressor matrix
      (original regressor matrix post-multiplied by \code{restrict.regMat})
      rather than the coefficients of the original regressor matrix
      are returned.}
   \item{\dots}{other arguments.}
}

\value{
   \code{coef.systemfit} returns a vector of all estimated coefficients.

   \code{coef.systemfit.equation} returns a vector of the estimated coefficients
      of a single equation.

   \code{coef.summary.systemfit} returns a matrix of all estimated coefficients,
      their standard errors, t-values, and p-values.

   \code{coef.summary.systemfit.equation} returns a matrix of the estimated
      coefficients of a single equation, their standard errors, t-values,
      and p-values.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{
   \code{\link{systemfit}}, \code{\link{coef}}
}


\examples{
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
system <- list( demand = eqDemand, supply = eqSupply )

## perform OLS on each of the equations in the system
fitols <- systemfit( system, data = Kmenta )

## all coefficients
coef( fitols )
coef( summary( fitols ) )

## coefficients of the first equation
coef( fitols$eq[[1]] )
coef( summary( fitols$eq[[1]] ) )

## coefficients of the second equation
coef( fitols$eq[[2]] )
coef( summary( fitols$eq[[2]] ) )

## estimation with restriction by modifying the regressor matrix
modReg <- matrix( 0, 7, 6 )
colnames( modReg ) <- c( "demIntercept", "demPrice", "demIncome",
   "supIntercept", "supPrice2", "supTrend" )
modReg[ 1, "demIntercept" ] <- 1
modReg[ 2, "demPrice" ]     <- 1
modReg[ 3, "demIncome" ]    <- 1
modReg[ 4, "supIntercept" ] <- 1
modReg[ 5, "supPrice2" ]    <- 1
modReg[ 6, "supPrice2" ]    <- 1
modReg[ 7, "supTrend" ]     <- 1
fitols3 <- systemfit( system, data = Kmenta, restrict.regMat = modReg )
coef( fitols3, modified.regMat  = TRUE )
coef( fitols3 )
}

\keyword{models}



