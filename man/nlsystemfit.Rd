
%  $Id: nlsystemfit.Rd,v 1.5 2005/03/18 15:25:32 henningsena Exp $


\name{nlsystemfit}
\alias{nlsystemfit}

\title{Nonlinear Equation System Estimation}

\description{
  Fits a set of structural nonlinear equations using Ordinary Least
  Squares (OLS), Seemingly Unrelated Regression (SUR),
  Two-Stage Least Squares (2SLS), Three-Stage Least Squares (3SLS).
}


\usage{ nlsystemfit( method="OLS", eqns, startvals,
                     eqnlabels=c(as.character(1:length(eqns))), inst=NULL,
                     data=list(), solvtol=.Machine$double.eps,
                     maxiter=1000, ... )
}


\arguments{
  \item{method}{the estimation method, one of "OLS", "SUR", "2SLS",
    "3SLS".}
  \item{eqns}{a list of structural equations to be estimated.}
  \item{startvals}{a list of starting parameter values for the minimization.}
  \item{eqnlabels}{an optional list of character vectors of names
    for the equation labels.}
  \item{inst}{one-sided model formula specifying instrumental variables
    or a list of one-sided model formulas if different instruments should
    be used for the different equations (only needed for 2SLS, 3SLS and
    GMM estimations).}
  \item{data}{an optional data frame containing the variables in the model.
    By default the variables are taken from the environment from which
    nlsystemfit is called.}
  \item{solvtol}{tolerance for detecting linear dependencies in the columns
    of X in the \code{\link{qr}} function calls.}
  \item{maxiter}{the maximum number of iterations for the \code{\link{nlm}} function.}
  \item{...}{arguments passed to \code{\link{nlm}}.}
}

\details{

  The nlsystemfit function relies on \code{\link{nlm}} to perform the
  minimization of the objective functions and the \code{\link{qr}} set
  of functions.

  A system of nonlinear equations can be written as:

  \deqn{\epsilon_{t} = q( y_t, x_t, \theta )}
  \deqn{z_{t} = Z( x_t )}

  where \eqn{\epsilon_{t}} are the residuals from the y observations and
  the function evaluated at the parameter estimates.

  The objective functions for the methods are:

  \tabular{lccc}{
    %    \hline
    Method \tab Instruments \tab Objective Function \tab Covariance of
    \eqn{\theta}\cr %\hline
    OLS \tab no \tab \eqn{r'r} \tab \eqn{(X(diag(S)^{-1}\bigotimes
      I)X)^{-1}}\cr %\hline
    SUR \tab no \tab \eqn{r'(diag(S)_{OLS}^{-1}\bigotimes I)r} \tab
    \eqn{(X(S^{-1}\bigotimes I)X)^{-1}}\cr %\hline
    2SLS \tab yes \tab \eqn{r'(I \bigotimes W)r} \tab
    \eqn{(X(diag(S)^{-1}\bigotimes I)X)^{-1}}\cr %\hline
    3SLS \tab yes \tab \eqn{r'(S_{2SLS}^{-1} \bigotimes W)r} \tab
    \eqn{(X(diag(S)^{-1}\bigotimes W)X)^{-1}} %\hline
  }

  where, r is a column vector for the residuals for each equation, S is
  variance-covariance matrix between the equations
  (\eqn{\hat{\sigma}_{ij} = (\hat{e}_i' \hat{e}_j) /
    \sqrt{(T - k_i)*(T - k_j)}}), X is matrix of the
  partial derivates with respect to the parameters, W is a matrix of the
  instrument variables \eqn{Z(Z'Z)^{-1}Z}, Z is a matrix of the
  instrument variables, and I is an nxn identity matrix.

  The SUR and 3SLS methods requires two solutions. The first solution
  for the SUR is an OLS solution to obtain the variance-covariance
  matrix. The 3SLS uses the variance-covatiance from a 2SLS solution,
  then fits all the equations simultaneously.

  The user should be aware that the function is \bold{VERY} sensative to
  the starting values and the nlm function may not converge. The nlm
  function will be called with the \code{typsize} argument set the
  absolute values of the starting values for the OLS and 2SLS
  methods. For the SUR and 3SLS methods, the \code{typsize} argument is
  set to the absolute values of the resulting OLS and 2SLS parameter
  estimates from the nlm result structre. In addition, the starting
  values for the SUR and 3SLS methods are obtained from the OLS and 2SLS
  parameter estimates to shorten the number of iterations. The number of
  iterations reported in the summary are only those used in the last
  call to nlm, thus the number of iterations in the OLS portion of the
  SUR fit and the 2SLS portion of the 3SLS fit are not included.
}

\value{
  \code{nlsystemfit} returns a list of the class \code{nlsystemfit.system} and
  contains all results that belong to the whole system.
  This list contains one special object: "eq". It is a list and contains
  one object for each estimated equation. These objects are of the class
  \code{nlsystemfit.equation} and contain the results that belong only to the
  regarding equation.

  The objects of the class \code{nlsystemfit.system} and
  \code{nlsystemfit.equation} have the following components (the elements of
  the latter are marked with an asterisk (\eqn{*})):

  \item{eq}{a list object that contains a list object for each equation.}
  \item{method}{estimation method.}
  \item{resids}{an \eqn{n \times g} matrix of the residuals.}
  \item{g}{number of equations.}
  \item{n}{total number of observations.}
  \item{k}{total number of coefficients.}
  \item{b}{vector of all estimated coefficients.}
  \item{se}{estimated standard errors of \code{b}.}
  \item{t}{t values for \code{b}.}
  \item{p}{p values for \code{b}.}
  \item{bcov}{estimated covariance matrix of \code{b}.}
  \item{rcov}{estimated residual covariance matrix.}
  \item{drcov}{determinant of \code{rcov}.}
  \item{rcovest}{residual covariance matrix used for estimation (only SUR and 3SLS).}
  \item{rcor}{estimated residual correlation matrix.}
  \item{nlmest}{results from the nlm function call}
  \item{solvetol}{tolerance level when inverting a matrix or calculating a determinant.}

  ## elements of the class nlsystemfit.eq
  \item{eq}{a list that contains the results that belong to the individual equations.}
  \item{eqnlabel*}{the equation label of the ith equation (from the labels list).}
  \item{formula*}{model formula of the ith equation.}
  \item{n*}{number of observations of the ith equation.}
  \item{k*}{number of coefficients/regressors in the ith equation.}
  \item{df*}{degrees of freedom of the ith equation.}
  \item{b*}{estimated coefficients of the ith equation.}
  \item{se*}{estimated standard errors of \code{b}.}
  \item{t*}{t values for \code{b}.}
  \item{p*}{p values for \code{b}.}
  \item{covb*}{estimated covariance matrix of \code{b}.}
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

  Gallant, R. H. (1987)
  \emph{Nonlinear Equation Estimation}, John Wiley and Sons, 610 pp.

  SAS Institute (1999)
  \emph{SAS/ETS User's Guide, Version 8}, Cary NC: SAS Institute 1546 pp.

}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com} }

\seealso{\code{\link{systemfit}}, \code{\link{nlm}}, and \code{\link{qr}}}

\examples{
library( systemfit )
data( ppine )

hg.formula <- hg ~ exp( h0 + h1*log(tht) + h2*tht^2 + h3*elev + h4*cr)
dg.formula <- dg ~ exp( d0 + d1*log(dbh) + d2*hg + d3*cr + d4*ba  )
labels <- list( "height.growth", "diameter.growth" )
inst <- ~ tht + dbh + elev + cr + ba
start.values <- c(h0=-0.5, h1=0.5, h2=-0.001, h3=0.0001, h4=0.08,
                  d0=-0.5, d1=0.009, d2=0.25, d3=0.005, d4=-0.02 )
model <- list( hg.formula, dg.formula )

model.ols <- nlsystemfit( "OLS", model, start.values, data=ppine, eqnlabels=labels )
print( model.ols )

model.sur <- nlsystemfit( "SUR", model, start.values, data=ppine, eqnlabels=labels )
print( model.sur )

model.2sls <- nlsystemfit( "2SLS", model, start.values, data=ppine,
   eqnlabels=labels, inst=inst )
print( model.2sls )

model.3sls <- nlsystemfit( "3SLS", model, start.values, data=ppine,
                                    eqnlabels=labels, inst=inst )
print( model.3sls )
}

\keyword{models}
\keyword{regression}
\keyword{nonlinear}



