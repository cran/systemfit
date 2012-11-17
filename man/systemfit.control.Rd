\name{systemfit.control}
\alias{systemfit.control}

\title{Create list of control parameters for systemfit}

\description{
   Create a list of control pararameters for function \code{\link{systemfit}}.
   All control parameters that are not passed to this function
   are set to default values.
}

\usage{
systemfit.control(
      maxiter = 1,
      tol = 1e-5,
      methodResidCov = "geomean",
      centerResiduals = FALSE,
      residCovRestricted = TRUE,
      residCovWeighted = FALSE,
      method3sls = "GLS",
      singleEqSigma = NULL,
      useMatrix = TRUE,
      solvetol = .Machine$double.eps,
      model = TRUE,
      x = FALSE,
      y = FALSE,
      z = FALSE )
}

\arguments{
   \item{maxiter}{maximum number of iterations for WLS, SUR, W2SLS and
      3SLS estimations.}
   \item{tol}{tolerance level indicating when to stop the iteration (only
      WLS, SUR, W2SLS and 3SLS estimations).}
   \item{methodResidCov}{method for calculating the estimated residual covariance
      matrix, one of "noDfCor", "geomean", "max", or "Theil" (see details).}
   \item{centerResiduals}{logical. Subtract the means from the residuals
      of each equation before calculating the estimated residual covariance matrix.}
   \item{residCovRestricted}{logical. If 'FALSE'
      the residual covariance matrix for a WLS, SUR, W2SLS, or 3SLS estimation
      is obtained from an unrestricted first-step estimation.}
   \item{residCovWeighted}{logical. If 'TRUE'
      the residual covariance matrix for a SUR or 3SLS estimation
      is obtained from a WLS or W2SLS estimation.}
   \item{method3sls}{method for calculating the 3SLS estimator,
      one of "GLS", "IV", "GMM", "Schmidt", or "EViews" (see details).}
   \item{singleEqSigma}{logical. use different \eqn{\sigma^2}s for each
      single equation to calculate the covariance matrix and the
      standard errors of the coefficients (only OLS and 2SLS)?
      If \code{singleEqSigma} is \code{NULL}, it is automatically determined:
      It is set to \code{TRUE}, if restrictions on the coefficients are imposed,
      and it is set to \code{FALSE} otherwise.}
   \item{useMatrix}{logical. Use package \code{Matrix} for matrix calculations?}
   \item{solvetol}{tolerance level for detecting linear dependencies
      when inverting a matrix or calculating a determinant (see
      \code{\link{solve}} and \code{\link{det}}).}
   \item{model, x, y, z}{
      logical. If 'TRUE' the corresponding components of the fit
      (the model frame, the model matrix, the response, and
      the matrix of instruments, respectively) are returned.}
}

\details{
   If the estimation is iterated
      (WLS, SUR, W2SLS or 3SLS estimation with \code{maxiter}>1),
      the convergence criterion is
   \deqn{\sqrt{ \frac{ \sum_i (b_{i,g} - b_{i,g-1})^2 }{ \sum_i b_{i,g-1}^2 }}
      < \code{tol}}
   (\eqn{b_{i,g}} is the ith coefficient of the gth iteration step).

   The method for calculating the estimated covariance matrix of the residuals
   (\eqn{\hat{\Sigma}}) can be one of the following
   (see Judge et al., 1985, p. 469): \cr
   if methodResidCov='noDfCor':
      \deqn{\hat{\sigma}_{ij} = \frac{\hat{e}_i' \hat{e}_j}{T}}
   if methodResidCov='geomean':
      \deqn{\hat{\sigma}_{ij} = \frac{\hat{e}_i' \hat{e}_j}
         {\sqrt{(T - k_i)*(T - k_j)}}}
   if methodResidCov='Theil':
      \deqn{\hat{\sigma}_{ij} = \frac{\hat{e}_i' \hat{e}_j}{T - k_i - k_j
      + tr[X_i(X_i'X_i)^{-1}X_i'X_j(X_j'X_j)^{-1}X_j']}}
   if methodResidCov='max':
      \deqn{\hat{\sigma}_{ij} = \frac{\hat{e}_i' \hat{e}_j}
         {T - \max( k_i, k_j)}}
   If \eqn{ i = j}, the formulas 'geomean', 'Theil', and 'max' are equal.
   All these three formulas yield unbiased estimators
   for the diagonal elements of the residual covariance matrix.
   If \eqn{i \neq j}, only formula 'Theil' yields an unbiased estimator for the residual
   covariance matrix, but it is not neccessarily positive semidefinit. Thus, it is
   doubtful whether formula 'Theil' is really superior to formula 'noDfCor'
   (Theil, 1971, p. 322).

   The methods for calculating the 3SLS estimator lead to identical results
   if the same instruments are used in all equations. If different instruments
   are used in the different equations, only the GMM-3SLS estimator ("GMM")
   and the 3SLS estimator proposed by Schmidt (1990) ("Schmidt") are consistent,
   whereas "GMM" is efficient relative to "Schmidt" (see Schmidt, 1990).

   If \code{residCovWeighted} is \code{TRUE},
   \code{\link{systemfit}} does a OLS or 2SLS estimation in a first step.
   It uses the residuals from the first-step estimation
   to calculate the residual covariance matrix
   that is used in a second-step WLS or W2SLS estimation.
   Then, it uses the residuals from the second-step estimation
   to calculate the residual covariance matrix
   that is used in a final SUR or 3SLS estimation.
   This three-step method is the default method of command "TSCS"
   in the software LIMDEP that carries out "SUR" estimations
   in which all coefficient vectors are constrained to be equal
   (personal information from W.H. Greene, 2006/02/16).
   If no cross-equation restrictions are imposed,
   \code{residCovWeighted} has no effect on the estimation results.
}

\value{
   A list of the above components.
}

\references{
  Judge, George G.; W. E. Griffiths; R. Carter Hill; Helmut Luetkepohl and Tsoung-Chao Lee (1985)
  \emph{The Theory and Practice of Econometrics, Second Edition}, Wiley.

  Schmidt, P. (1990)
  \emph{Three-Stage Least Squares with different Instruments for different equations},
  Journal of Econometrics 43, p. 389-394.

  Theil, H. (1971)
  \emph{Principles of Econometrics}, Wiley, New York.
}

\author{
  Arne Henningsen \email{arne.henningsen@googlemail.com}
}

\seealso{\code{\link{systemfit}}}

\examples{
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
eqSystem <- list( demand = eqDemand, supply = eqSupply )

## SUR estimation: calculation of residual covariance
## matrix without correction for degrees of freedom
fitsur <- systemfit( eqSystem, "SUR", data = Kmenta,
    control = systemfit.control( methodResidCov = "noDfCor" ) )
print( fitsur )
}

\keyword{models}
\keyword{regression}



