
%    $Id: systemfit.Rd 1031 2009-06-17 16:14:50Z arne $


\name{systemfit}
\alias{systemfit}

\title{Linear Equation System Estimation}

\description{
  Fits a set of linear structural equations using Ordinary Least
  Squares (OLS), Weighted Least Squares (WLS), Seemingly Unrelated Regression (SUR),
  Two-Stage Least Squares (2SLS), Weighted Two-Stage Least Squares (W2SLS)
  or Three-Stage Least Squares (3SLS).
}

\usage{
systemfit( formula, method = "OLS",
           inst=NULL, data=list(),
           restrict.matrix = NULL, restrict.rhs = NULL, restrict.regMat = NULL,
           pooled = FALSE, control = systemfit.control( ... ), ... )
}

\arguments{
   \item{formula}{an object of class \code{formula}
      (for single-equation models)
      or (typically) a list of objects of class \code{formula}
      (for multiple-equation models);
      if argument \code{data} is of class \code{plm.dim}
      (created with \code{plm.data}),
      this argument must be a single object of class \code{formula}
      that represents the formula to be estimated for all individuals.}
   \item{method}{the estimation method, one of "OLS", "WLS", "SUR",
      "2SLS", "W2SLS", or "3SLS" (see details);
      iterated estimation methods can be specified by setting control parameter
      \code{maxiter} larger than 1 (e.g. 500).}
   \item{inst}{one-sided model formula specifying instrumental variables
      or a list of one-sided model formulas if different instruments should
      be used for the different equations (only needed for 2SLS, W2SLS, and
      3SLS estimations).}
   \item{data}{an optional data frame
      containing the variables in the model.
      By default the variables are taken from the environment from which
      systemfit is called.}
   \item{restrict.matrix}{an optional j x k matrix to impose linear
      restrictions on the coefficients by \code{restrict.matrix} * \eqn{b} = \code{restrict.rhs}
      (j = number of restrictions, k = number of all coefficients,
      \eqn{b} = vector of all coefficients)
      or a character vector giving the restrictions in symbolic form
      (see documentation of \code{\link[car]{linear.hypothesis}}
      in package "car" for details).
      The number and the names of the coefficients can be obtained
      by estimating the system without restrictions and applying
      the \code{coef} method to the returned object.}
   \item{restrict.rhs}{an optional vector with j elements to impose linear
      restrictions (see \code{restrict.matrix}); default is a vector
      that contains j zeros.}
   \item{restrict.regMat}{an optional matrix to impose restrictions
      on the coefficients by post-multiplying the regressor matrix
      with this matrix (see details).}
   \item{control}{list of control parameters.
      The default is constructed by the function \code{\link{systemfit.control}}.
      See the documentation of \code{\link{systemfit.control}} for details.}
   \item{pooled}{logical, restrict coefficients to be equal in all equations
      (only for panel-like data).}
   \item{...}{arguments passed to \code{\link{systemfit.control}}.}
}

\details{
   The estimation of systems of equations with unequal numbers
   of observations has not been thoroughly tested yet.
   Currently, \code{systemfit} calculates the residual covariance matrix
   only from the residuals/observations that are available
   in all equations.

   If argument \code{data} is of class \code{plm.dim}
   (created with \code{plm.data} and thus, contains panel data in long format),
   argument \code{formula} must be a single equation
   that is applied to all individuals.
   In this case, argument \code{pooled} specifies
   whether the coefficients are restricted to be equal for all
   individuals.

   If argument \code{restrict.regMat} is specified,
   the regressor matrix \eqn{X} is post-multiplied by this matrix:
   \eqn{X^{*} = X \cdot} \code{restrict.regMat}.
   Then, this modified regressor matrix \eqn{X^{*}} is used for the
   estimation of the coefficient vector \eqn{b^{*}}.
   This means that the coefficients of the original regressors (\eqn{X}),
   vector \eqn{b},
   can be represented by \eqn{b =} \code{restrict.regMat} \eqn{\cdot b^{*}}.
   If \code{restrict.regMat} is a non-singular quadratic matrix,
   there are no restrictions on the coefficients imposed,
   but the coefficients \eqn{b^{*}} are linear combinations of
   the original coefficients \eqn{b}.
   If \code{restrict.regMat} has less columns than rows,
   linear restrictions are imposed on the coefficients \eqn{b}.
   However, imposing linear restrictions
   by the \code{restrict.regMat} matrix is less flexible than
   by providing the matrix \code{restrict.matrix} and the
   vector \code{restrict.rhs}.
   The advantage of imposing restrictions on the coefficients
   by the matrix \code{restrict.regMat} is that the matrix,
   which has to be inverted during the estimation,
   gets smaller by this procedure, while it gets larger
   if the restrictions are imposed by \code{restrict.matrix} and \code{restrict.rhs}.

   In the context of multi-equation models, 
   the term \dQuote{weighted} in \dQuote{weighted least squares} (WLS) and
   \dQuote{weighted two-stage least squares} (W2SLS) means
   that the \emph{equations} might have different weights
   and \emph{not}
   that the \emph{observations} have different weights.
}

\value{
  \code{systemfit} returns a list of the class \code{systemfit} and
  contains all results that belong to the whole system.
  This list contains one special object: "eq". It is a list and contains
  one object for each estimated equation. These objects are of the class
  \code{systemfit.equation} and contain the results that belong only to the
  regarding equation.

  The objects of the class \code{systemfit} and
  \code{systemfit.equation} have the following components (the elements of
  the latter are marked with an asterisk (\eqn{*})):

  \item{call}{the matched call.}
  \item{method}{estimation method.}
  \item{rank}{total number of linear independent coefficients
      = number of coefficients minus number of linear restrictions.}
  \item{df.residual}{degrees of freedom of the whole system.}
  \item{iter}{number of iteration steps.}
  \item{coefficients}{vector of all estimated coefficients.}
  \item{coefCov}{estimated covariance matrix of \code{coefficients}.}
  \item{residCov}{estimated residual covariance matrix.}
  \item{residCovEst}{residual covariance matrix used for estimation
      (only WLS, W2SLS, SUR and 3SLS).}
  \item{restrict.matrix}{the restriction matrix.}
  \item{restrict.rhs}{the restriction vector.}
  \item{restrict.regMat}{matrix used to impose restrictions
      on the coefficients by post-multiplying the regressor matrix
      with this matrix.}
  \item{control}{list of control parameters used for the estimation.}
  \item{panelLike}{logical. Was this an analysis with panel-like data?}

  ## elements of the class systemfit.eq
  \item{eq}{a list that contains the results that belong to the individual equations.}
  \item{eqnLabel*}{the label of this equation.}
  \item{eqnNo*}{the number of this equation.}
  \item{terms*}{the 'terms' object used for the ith equation.}
  \item{inst*}{instruments of the ith equation (only 2SLS, W2SLS, and 3SLS).}
  \item{rank*}{number of linear independent coefficients in the ith equation
     (differs from the number of coefficients only if there are
     restrictions that are not cross-equation).}
  \item{nCoef.sys*}{total number of coefficients in all equations.}
  \item{rank.sys*}{total number of linear independent coefficients
      in all equations.}
  \item{df.residual*}{degrees of freedom of the ith equation.}
  \item{df.residual.sys*}{degrees of freedom of the whole system.}
  \item{coefficients*}{estimated coefficients of the ith equation.}
  \item{covb*}{estimated covariance matrix of \code{coefficients}.}
  \item{model*}{if requested (the default), the model frame of the ith equation.}
  \item{x*}{if requested, the model matrix of the ith equation.}
  \item{y*}{if requested, the response of the ith equation.}
  \item{z*}{if requested, the matrix of instrumental variables
      of the ith equation (only 2SLS, W2SLS, and 3SLS).}
  \item{fitted.values*}{vector of fitted values of the ith equation.}
  \item{residuals*}{vector of residuals of the ith equation.}
}

\references{

  Greene, W. H. (2003)
  \emph{Econometric Analysis, Fifth Edition}, Prentice Hall.

  Judge, George G.; W. E. Griffiths; R. Carter Hill; Helmut Luetkepohl and Tsoung-Chao Lee (1985)
  \emph{The Theory and Practice of Econometrics, Second Edition}, Wiley.

  Kmenta, J. (1997)
  \emph{Elements of Econometrics, Second Edition}, University of
  Michigan Publishing.

  Schmidt, P. (1990)
  \emph{Three-Stage Least Squares with different Instruments for different equations},
  Journal of Econometrics 43, p. 389-394.

  Theil, H. (1971)
  \emph{Principles of Econometrics}, Wiley, New York.
}

\author{Arne Henningsen \email{arne.henningsen@googlemail.com},\cr
  Jeff D. Hamann \email{jeff.hamann@forestinformatics.com}
}

\seealso{\code{\link{lm}} and \code{\link{nlsystemfit}}}

\examples{
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
system <- list( demand = eqDemand, supply = eqSupply )

## OLS estimation
fitols <- systemfit( system, data=Kmenta )
print( fitols )

## OLS estimation with 2 restrictions
Rrestr <- matrix(0,2,7)
Rrestr[1,3] <-  1
Rrestr[1,7] <- -1
Rrestr[2,2] <- -1
Rrestr[2,5] <-  1
qrestr <- c( 0, 0.5 )
fitols2 <- systemfit( system, data = Kmenta,
                      restrict.matrix = Rrestr, restrict.rhs = qrestr )
print( fitols2 )

## OLS estimation with the same 2 restrictions in symbolic form
restrict <- c( "demand_income - supply_trend = 0",
   "- demand_price + supply_price = 0.5" )
fitols2b <- systemfit( system, data = Kmenta, restrict.matrix = restrict )
print( fitols2b )

# test whether both restricted estimators are identical
all.equal( coef( fitols2 ), coef( fitols2b ) )

## OLS with restrictions on the coefficients by modifying the regressor matrix
## with argument restrict.regMat
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
print( fitols3 )


## iterated SUR estimation
fitsur <- systemfit( system, "SUR", data = Kmenta, maxit = 100 )
print( fitsur )

## 2SLS estimation
inst <- ~ income + farmPrice + trend
fit2sls <- systemfit( system, "2SLS", inst = inst, data = Kmenta )
print( fit2sls )

## 2SLS estimation with different instruments in each equation
inst1 <- ~ income + farmPrice
inst2 <- ~ income + farmPrice + trend
instlist <- list( inst1, inst2 )
fit2sls2 <- systemfit( system, "2SLS", inst = instlist, data = Kmenta )
print( fit2sls2 )

## 3SLS estimation with GMM-3SLS formula
inst <- ~ income + farmPrice + trend
fit3sls <- systemfit( system, "3SLS", inst = inst, data = Kmenta,
   method3sls = "GMM" )
print( fit3sls )


## Examples how to use systemfit() with panel-like data
## Repeating the SUR estimations in Greene (2003, p. 351)
data( "GrunfeldGreene" )
library( plm )
GGPanel <- plm.data( GrunfeldGreene, c( "firm", "year" ) )
formulaGrunfeld <- invest ~ value + capital
# SUR
greeneSur <- systemfit( formulaGrunfeld, "SUR",
   data = GGPanel, methodResidCov = "noDfCor" )
summary( greeneSur )
# SUR Pooled
greeneSurPooled <- systemfit( formulaGrunfeld, "SUR",
   data = GGPanel, pooled = TRUE, methodResidCov = "noDfCor",
   residCovWeighted = TRUE )
summary( greeneSurPooled )

## Further examples are in the documentation to the data sets
## 'KleinI' and 'GrunfeldGreene'.
}

\keyword{models}
\keyword{regression}



