
% 	$Id: systemfit.Rd,v 1.7 2004/01/30 20:03:43 hamannj Exp $	


\name{systemfit}
\alias{systemfit}

\title{Equation System Estimation}

\description{
  Fits a set of linear structural equations using Ordinary Least
  Squares (OLS), Weighted Least Squares (WLS), Seemingly Unrelated Regression (SUR),
  Two-Stage Least Squares (2SLS), Weighted Two-Stage Least Squares (W2SLS)
  or Three-Stage Least Squares (3SLS).
}

\usage{
systemfit( method, eqns, eqnlabels=c(as.character(1:length(eqns))),
           inst=NULL, data=list(), R.restr=NULL,
           q.restr=matrix(0,max(nrow(R.restr),0),1),
           TX=NULL, maxiter=1, tol=1e-5,
           rcovformula=1, formula3sls="GLS",
           probdfsys=!(is.null(R.restr) & is.null(TX)),
           single.eq.sigma=(is.null(R.restr) & is.null(TX)),
           solvetol=.Machine$double.eps )
}

\arguments{
  \item{method}{the estimation method, one of "OLS", "WLS", "SUR",
    "2SLS", "W2SLS" or "3SLS".}
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
}

\details{
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
     the convergence criterion is \cr
   \eqn{\sqrt{\sum_i (b_{i,g} - b_{i,g-1})^2 \left/ \sum_i b_{i,g-1}^2 \right.}}
   < \code{tol} (\eqn{b_{i,g}} is the ith coefficient of the gth iteration step).

   The formula to calculate the estimated covariance matrix of the residuals
   (\eqn{\hat{\Sigma}}) can be one of the following (see Judge et al., 1985, p. 469): \cr
   if rcovformula=0: \eqn{\hat{\sigma}_{ij} = (\hat{e}_i' \hat{e}_j) / T}; \cr
   if rcovformula=1: \eqn{\hat{\sigma}_{ij} = (\hat{e}_i' \hat{e}_j) /
       \sqrt{(T - k_i)*(T - k_j)}}; \cr
   if rcovformula=2: \eqn{\hat{\sigma}_{ij} = (\hat{e}_i' \hat{e}_j) / (T - k_i - k_j
   + tr[(X_i'X_i)^{-1}X_i'X_j(X_j'X_j)^{-1}X_j'X_i]}. \cr
   If \eqn{k_i = k_j}, formula 1 and 2 are equal and yield an unbiased estimator for the
   residual covariance matrix.
   If \eqn{k_i \neq k_j}, only formula 2 yields an unbiased estimator for the residual
   covariance matrix, but it is not neccessarily positive semidefinit and
   its inverse is \bold{not} an unbiased estimator for the inverse of the residual
   covariance matrix. Thus, it is doubtful whether formula 2 is really superior
   to formula 1 (see Theil, 1971, p. 322).

   The formulas to calculate the 3SLS estimator lead to identical results
   if the same instruments are used in all equations. If different instruments
   are used in the different equations, only the GMM-3SLS estimator ("GMM")
   and the 3SLS estimator proposed by Schmidt (1990) ("Schmidt") are consistent,
   whereas "GMM" is efficient relative to "Schmidt" (see Schmidt, 1990).
}

\value{
  \code{systemfit} returns a list of the class \code{systemfit.system} and
  contains all results that belong to the whole system.
  This list contains one special object: "eq". It is a list and contains
  one object for each estimated equation. These objects are of the class
  \code{systemfit.equation} and contain the results that belong only to the
  regarding equation.

  The objects of the class \code{systemfit.system} and
  \code{systemfit.equation} have the following components (the elements of
  the latter are marked with an asterisk (\eqn{*})):

  \item{method}{estimation method.}
  \item{g}{number of equations.}
  \item{n}{total number of observations.}
  \item{k}{total number of coefficients.}
  \item{ki}{total number of linear independent coefficients.}
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
  \item{rcor}{estimated residual correlation matrix.}
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
  \item{predicted*}{vector of predicted values of the ith equation.}
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

  Greene, W. H. (1993)
  \emph{Econometric Analysis, Second Edition}, Macmillan.

  Greene, W. H. (2002)
  \emph{Econometric Analysis, Fifth Edition}, Prentice Hall.

  Judge, George G.; W. E. Griffiths; R. Carter Hill; Helmut Lütkepohl and Tsoung-Chao Lee (1985)
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

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com}
  Arne Henningsen (\url{http://www.arne-henningsen.de}) }
}

\seealso{\code{\link{lm}} and \code{\link{nlsystemfit}}}

\examples{
\dontrun{library( systemfit )}
data( kmenta )
demand <- q ~ p + d
supply <- q ~ p + f + a
labels <- list( "demand", "supply" )
system <- list( demand, supply )

## OLS estimation
fitols <- systemfit("OLS", system, labels, data=kmenta )
print( fitols )

## OLS estimation with 2 restrictions
Rrestr <- matrix(0,2,7)
qrestr <- matrix(0,2,1)
Rrestr[1,3] <-  1
Rrestr[1,7] <- -1
Rrestr[2,2] <- -1
Rrestr[2,5] <-  1
qrestr[2,1] <-  0.5
fitols2 <- systemfit("OLS", system, labels, data=kmenta,
                      R.restr=Rrestr, q.restr=qrestr )
print( fitols2 )

## iterated SUR estimation
fitsur <- systemfit("SUR", system, labels, data=kmenta, maxit=100 )
print( fitsur )

## 2SLS estimation
inst <- ~ d + f + a
fit2sls <- systemfit( "2SLS", system, labels, inst, kmenta )
print( fit2sls )

## 2SLS estimation with different instruments in each equation
inst1 <- ~ d + f
inst2 <- ~ d + f + a
instlist <- list( inst1, inst2 )
fit2sls2 <- systemfit( "2SLS", system, labels, instlist, kmenta )
print( fit2sls2 )

## 3SLS estimation with GMM-3SLS formula
inst <- ~ d + f + a
fit3sls <- systemfit( "3SLS", system, labels, inst, kmenta, formula3sls="GMM" )
print( fit3sls )

}

\keyword{models}
\keyword{regression}



