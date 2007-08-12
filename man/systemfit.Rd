
%    $Id: systemfit.Rd 479 2007-08-12 08:13:50Z henningsena $


\name{systemfit}
\alias{systemfit}
\alias{systemfitClassic}

\title{Linear Equation System Estimation}

\description{
  Fits a set of linear structural equations using Ordinary Least
  Squares (OLS), Weighted Least Squares (WLS), Seemingly Unrelated Regression (SUR),
  Two-Stage Least Squares (2SLS), Weighted Two-Stage Least Squares (W2SLS)
  or Three-Stage Least Squares (3SLS).
}

\usage{
systemfit( method, eqns, eqnlabels=names(eqns),
           inst=NULL, data=list(), R.restr=NULL,
           q.restr=matrix(0,max(nrow(R.restr),0),1),
           TX=NULL, maxiter=1, tol=1e-5,
           rcovformula=1, centerResiduals = FALSE, formula3sls="GLS",
           probdfsys=!(is.null(R.restr) & is.null(TX)),
           single.eq.sigma=(is.null(R.restr) & is.null(TX)),
           solvetol=.Machine$double.eps,
           saveMemory=( nrow(data) * length(eqns) > 1000 &&
              length(data) > 0 ) )

systemfitClassic( method, formula, eqnVar, timeVar, data,
                  pooled = FALSE, ... )
}

\arguments{
   \item{method}{the estimation method, one of "OLS", "WLS", "SUR", "WSUR",
      "2SLS", "W2SLS", "3SLS", or "W3SLS" (see details);
      iterated estimation methods can be specified by setting argument
      \code{maxiter} larger than 1 (e.g. 500).}
   \item{eqns}{a list of structural equations to be estimated; a regression
      constant is implied if not explicitly omitted.}
   \item{eqnlabels}{an optional list of character vectors of names
      for the equation labels.}
   \item{inst}{one-sided model formula specifying instrumental variables
      or a list of one-sided model formulas if different instruments should
      be used for the different equations (only needed for 2SLS, W2SLS and
      3SLS estimations).}
   \item{data}{an optional data frame containing the variables in the model.
      By default the variables are taken from the environment from which
      systemfit is called.}
   \item{R.restr}{an optional j x k matrix to impose linear
      restrictions on the parameters by \code{R.restr} * \eqn{b} = \code{q.restr}
      (j = number of restrictions, k = number of all parameters,
      \eqn{b} = vector of all parameters).}
   \item{q.restr}{an optional j x 1 matrix to impose linear
      restrictions (see \code{R.restr}); default is a j x 1 matrix
      that contains only zeros.}
   \item{TX}{an optional matrix to transform the regressor matrix and,
      hence, also the coefficient vector (see details).}
   \item{maxiter}{maximum number of iterations for WLS, SUR, W2SLS and
      3SLS estimations.}
   \item{tol}{tolerance level indicating when to stop the iteration (only
      WLS, SUR, W2SLS and 3SLS estimations).}
   \item{rcovformula}{formula to calculate the estimated residual covariance
      matrix (see details).}
   \item{centerResiduals}{logical. Subtract the means from the residuals
      of each equation before calculating the estimated residual covariance matrix.}
   \item{formula3sls}{formula for calculating the 3SLS estimator,
      one of "GLS", "IV", "GMM", "Schmidt" or "EViews" (see details).}
   \item{probdfsys}{use the degrees of freedom of the whole system
      (in place of the degrees of freedom of the single equation)
      to calculate prob values for the t-test of individual parameters.}
   \item{single.eq.sigma}{use different \eqn{\sigma^2}s for each
      single equation to calculate the covariance matrix and the
      standard errors of the coefficients (only OLS and 2SLS).}
   \item{solvetol}{tolerance level for detecting linear dependencies
      when inverting a matrix or calculating a determinant (see
      \code{\link{solve}} and \code{\link{det}}).}
   \item{saveMemory}{logical. Save memory by omitting some calculation that
      are not crucial for the basic estimation (e.g McElroy's
      \eqn{R^2})?}
   \item{formula}{formula to be estimated (for each equation).}
   \item{eqnVar}{variable name indicating the equation to which the
      observation belongs.}
   \item{timeVar}{variable name indicating the time.}
   \item{pooled}{logical, restrict coefficients to be equal in all equations.}
   \item{...}{arguments passed to \code{\link{systemfit}}.}
}

\details{
   \code{systemfitClassic} is a wrapper function for \code{systemfit}
   that can be applied to panel-like data in long format
   if the regressors are the same for all equations.

   If argument \code{method} is "WSUR" or "W3SLS",
   the "SUR" or "3SLS" estimation uses a residual variance covariance matrix
   that is calculated from a "WLS" or "W2SLS" estimation, respectively
   (and not from an "OLS" or "2SLS" estimation as for a standard "SUR" or "3SLS"
   estimation).
   The "WSUR" method is the default method of command "TSCS"
   in the software LIMDEP that carries out "SUR" estimations
   in which all coefficient vectors are constrained to be equal
   (personal information from W.H. Greene, 2006/02/16).
   If no cross-equation restrictions are imposed, "WSUR" and "W3SLS"
   generate identical results compared to "SUR" and "3SLS", respectively.

   The matrix \code{TX} transforms the regressor matrix (\eqn{X}) by
   \eqn{X^{*} = X *} \code{TX}. Thus, the vector of coefficients is now
   \eqn{b =} \code{TX} \eqn{\cdot b^{*}} , where \eqn{b} is the original (stacked) vector
   of all coefficients and \eqn{b^{*}} is the new coefficient vector that is
   estimated instead. Thus, the elements of vector \eqn{b} are
   \eqn{b_i = \sum_j TX_{ij} \cdot b^{*}_j} \cr
   The \code{TX} matrix can be used to change the order of the
   coefficients and also to restrict coefficients (if \code{TX} has less
   columns than it has rows). However restricting coefficients
   by the \code{TX} matrix is less powerfull and flexible than the
   restriction by providing the \code{R.restr} matrix and the
   \code{q.restr} vector. The advantage of restricting the coefficients
   by the \code{TX} matrix is that the matrix that is inverted for
   estimation gets smaller by this procedure, while it gets larger
   if the restrictions are imposed by \code{R.restr} and \code{q.restr}.

   If iterated (WLS, SUR, W2SLS or 3SLS estimation with \code{maxit}>1),
      the convergence criterion is
   \deqn{\sqrt{ \frac{ \sum_i (b_{i,g} - b_{i,g-1})^2 }{ \sum_i b_{i,g-1}^2 }}
      < \code{tol}}
   (\eqn{b_{i,g}} is the ith coefficient of the gth iteration step).

   The formula to calculate the estimated covariance matrix of the residuals
   (\eqn{\hat{\Sigma}}) can be one of the following
   (see Judge et al., 1985, p. 469): \cr
   if rcovformula=0:
      \deqn{\hat{\sigma}_{ij} = \frac{\hat{e}_i' \hat{e}_j}{T}}
   if rcovformula=1 or rcovformula='geomean':
      \deqn{\hat{\sigma}_{ij} = \frac{\hat{e}_i' \hat{e}_j}
         {\sqrt{(T - k_i)*(T - k_j)}}}
   if rcovformula=2 or rcovformula='Theil':
      \deqn{\hat{\sigma}_{ij} = \frac{\hat{e}_i' \hat{e}_j}{T - k_i - k_j
      + tr[X_i(X_i'X_i)^{-1}X_i'X_j(X_j'X_j)^{-1}X_j']}}
   if rcovformula=3 or rcovformula='max':
      \deqn{\hat{\sigma}_{ij} = \frac{\hat{e}_i' \hat{e}_j}
         {T - \max( k_i, k_j)}}
   If \eqn{i = j}, formula 1, 2 and 3 are equal.
   All these three formulas yield unbiased estimators
   for the diagonal elements of the residual covariance matrix.
   If \eqn{i \neq j}, only formula 2 yields an unbiased estimator for the residual
   covariance matrix, but it is not neccessarily positive semidefinit. Thus, it is
   doubtful whether formula 2 is really superior to formula 1
   (Theil, 1971, p. 322).

   The formulas to calculate the 3SLS estimator lead to identical results
   if the same instruments are used in all equations. If different instruments
   are used in the different equations, only the GMM-3SLS estimator ("GMM")
   and the 3SLS estimator proposed by Schmidt (1990) ("Schmidt") are consistent,
   whereas "GMM" is efficient relative to "Schmidt" (see Schmidt, 1990).
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

  \item{method}{estimation method.}
  \item{g}{number of equations.}
  \item{n}{total number of observations.}
  \item{k}{total number of coefficients.}
  \item{ki}{total number of linear independent coefficients.}
  \item{df}{degrees of freedom of the whole system.}
  \item{iter}{number of iteration steps.}
  \item{b}{vector of all estimated coefficients.}
  \item{bt}{coefficient vector transformed by \code{TX}.}
  \item{se}{estimated standard errors of \code{b}.}
  \item{t}{t values for \code{b}.}
  \item{p}{p values for \code{b}.}
  \item{bcov}{estimated covariance matrix of \code{b}.}
  \item{btcov}{covariance matrix of \code{bt}.}
  \item{rcov}{estimated residual covariance matrix.}
  \item{drcov}{determinant of \code{rcov}.}
  \item{rcovest}{residual covariance matrix used for estimation (only SUR and 3SLS).}
  \item{olsr2}{System OLS R-squared value.}
  \item{mcelr2}{McElroys R-squared value for the system (only SUR and 3SLS).}
  \item{y}{vector of all (stacked) endogenous variables}
  \item{x}{matrix of all (diagonally stacked) regressors}
  \item{h}{matrix of all (diagonally stacked) instrumental variables (only 2SLS and 3SLS)}
  \item{data}{data frame of the whole system (including instruments)}
  \item{R.restr}{the restriction matrix.}
  \item{q.restr}{the restriction vector.}
  \item{TX}{matrix used to transform the regressor matrix.}
  \item{maxiter}{maximum number of iterations.}
  \item{tol}{tolerance level indicating when to stop the iteration}
  \item{rcovformula}{formula to calculate the estimated residual covariance
    matrix}
  \item{formula3sls}{formula for calculating the 3SLS estimator.}
  \item{probdfsys}{system degrees of freedom to calculate prob values?.}
  \item{single.eq.sigma}{different \eqn{\sigma^2}s for each single equation?.}
  \item{solvetol}{tolerance level when inverting a matrix or calculating a determinant.}
  \item{data.name}{name of the data.frame used for estimation.}

  ## elements of the class systemfit.eq
  \item{eq}{a list that contains the results that belong to the individual equations.}
  \item{eqnlabel*}{the equation label of the ith equation (from the labels list).}
  \item{formula*}{model formula of the ith equation.}
  \item{inst*}{instruments of the ith equation (only 2SLS and 3SLS).}
  \item{n*}{number of observations of the ith equation.}
  \item{k*}{number of coefficients/regressors in the ith equation (including the constant).}
  \item{ki*}{number of linear independent coefficients in the ith equation (including the constant
              differs from \code{k} only if there are restrictions that are not cross-equation).}
  \item{df*}{degrees of freedom of the ith equation.}

  \item{b*}{estimated coefficients of the ith equation.}
  \item{se*}{estimated standard errors of \code{b}.}
  \item{t*}{t values for \code{b}.}
  \item{p*}{p values for \code{b}.}
  \item{covb*}{estimated covariance matrix of \code{b}.}

  \item{y*}{vector of endogenous variable (response values) of the ith equation.}
  \item{x*}{matrix of regressors (model matrix) of the ith equation.}
  \item{h*}{matrix of instrumental variables of the ith equation (only 2SLS and 3SLS).}
  \item{data*}{data frame (including instruments) of the ith equation.}
  \item{fitted*}{vector of fitted values of the ith equation.}
  \item{residuals*}{vector of residuals of the ith equation.}

  \item{ssr*}{sum of squared residuals of the ith equation.}
  \item{mse*}{estimated variance of the residuals (mean of squared errors) of the ith equation.}
  \item{s2*}{estimated variance of the residuals (\eqn{\hat{\sigma}^2}) of the ith equation.}
  \item{rmse*}{estimated standard error of the residulas (square root of mse) of the ith equation.}
  \item{s*}{estimated standard error of the residuals (\eqn{\hat{\sigma}}) of the ith equation.}

  \item{r2*}{R-squared (coefficient of determination).}
  \item{adjr2*}{adjusted R-squared value.}
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

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com},\cr
  Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de} }
}

\seealso{\code{\link{lm}} and \code{\link{nlsystemfit}}}

\examples{
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
system <- list( demand = eqDemand, supply = eqSupply )

## OLS estimation
fitols <- systemfit("OLS", system, data=Kmenta )
print( fitols )

## OLS estimation with 2 restrictions
Rrestr <- matrix(0,2,7)
qrestr <- matrix(0,2,1)
Rrestr[1,3] <-  1
Rrestr[1,7] <- -1
Rrestr[2,2] <- -1
Rrestr[2,5] <-  1
qrestr[2,1] <-  0.5
fitols2 <- systemfit("OLS", system, data = Kmenta,
                      R.restr = Rrestr, q.restr = qrestr )
print( fitols2 )

## iterated SUR estimation
fitsur <- systemfit("SUR", system, data = Kmenta, maxit = 100 )
print( fitsur )

## 2SLS estimation
inst <- ~ income + farmPrice + trend
fit2sls <- systemfit( "2SLS", system, inst = inst, data = Kmenta )
print( fit2sls )

## 2SLS estimation with different instruments in each equation
inst1 <- ~ income + farmPrice
inst2 <- ~ income + farmPrice + trend
instlist <- list( inst1, inst2 )
fit2sls2 <- systemfit( "2SLS", system, inst = instlist, data = Kmenta )
print( fit2sls2 )

## 3SLS estimation with GMM-3SLS formula
inst <- ~ income + farmPrice + trend
fit3sls <- systemfit( "3SLS", system, inst = inst, data = Kmenta,
   formula3sls = "GMM" )
print( fit3sls )


## Examples how to use systemfitClassic()
## Repeating the OLS and SUR estimations in Theil (1971, pp. 295, 300)
data( "GrunfeldTheil" )
formulaGrunfeld <- invest ~ value + capital
# OLS
theilOls <- systemfitClassic( "OLS", formulaGrunfeld, "firm", "year",
   data = GrunfeldTheil )
summary( theilOls )
# SUR
theilSur <- systemfitClassic( "SUR", formulaGrunfeld, "firm", "year",
   data = GrunfeldTheil, rcovformula = 0 )
summary( theilSur )


## Further examples are in the documentation to the data sets
## 'KleinI' and 'GrunfeldGreene'.
}

\keyword{models}
\keyword{regression}



