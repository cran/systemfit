\name{summary.systemfit}
\alias{summary.systemfit}
\alias{summary.systemfit.equation}
\alias{print.summary.systemfit}
\alias{print.summary.systemfit.equation}
\title{Summary of systemfit estimation}

\description{
   These functions create and print summary results
   of the estimated equation system.
}

\usage{
   \method{summary}{systemfit}( object, useDfSys = NULL,
      residCov = TRUE, equations = TRUE, ... )

   \method{summary}{systemfit.equation}( object, useDfSys = NULL, ... )

   \method{print}{summary.systemfit}( x,
      digits = max( 3, getOption("digits") - 1 ),
      residCov = x$printResidCov, equations = x$printEquations, ... )

   \method{print}{summary.systemfit.equation}( x,
      digits = max( 3, getOption("digits") - 1 ), ... )
}

\arguments{
   \item{object}{an object of class \code{systemfit} or
      \code{systemfit.equation}.}
   \item{x}{an object of class \code{summary.systemfit} or
      \code{summary.systemfit.equation}.}
   \item{useDfSys}{logical. Use the degrees of freedom of the whole system
      (in place of the degrees of freedom of the single equation)
      to calculate prob values for the t-test of individual coefficients.
      If it not specified (\code{NULL}), it is set to \code{TRUE}
      if restrictions on the coefficients are imposed and \code{FALSE} otherwise.}
   \item{digits}{number of digits to print.}
   \item{residCov}{logical. If \code{TRUE},
      the residual correlation matrix, the residual covariance matrix,
      and its determinant are printed.}
   \item{equations}{logical. If \code{TRUE},
      summary results of each equation are printed.
      If \code{FALSE}, just the coefficients are printed.}
   \item{...}{not used by user.}
}

\value{
   Applying \code{summary} on an object of class \code{systemfit}
   returns a list of class \code{summary.}\hspace{0pt}\code{systemfit}.
   Applying \code{summary} on an object of class
   \code{systemfit.}\hspace{0pt}\code{equation}
   returns a list of class
   \code{summary.}\hspace{0pt}\code{systemfit.}\hspace{0pt}\code{equation}.
   An object of class \code{summary.}\hspace{0pt}\code{systemfit}
   contains all results that belong to the whole system.
   This list contains one special object: \code{eq}.
   This is a list and contains objects of class
   \code{summary.}\hspace{0pt}\code{systemfit.}\hspace{0pt}\code{equation}.
   These objects contain the results that belong to each of the eatimated equations.

   The objects of classes \code{summary.}\hspace{0pt}\code{systemfit} and
   \code{summary.}\hspace{0pt}\code{systemfit.}\hspace{0pt}\code{equation}
   have the following components
   (elements that are marked with a \eqn{*} are available only in objects of
   class \code{summary.}\hspace{0pt}\code{systemfit};
   elements that are marked with a \eqn{+} are available only in objects of
   class \code{summary.}\hspace{0pt}\code{systemfit.}\hspace{0pt}\code{equation}):

   \item{method}{estimation method.}
   \item{residuals}{residuals.}
   \item{coefficients}{a matrix with columns for the estimated coefficients,
      their standard errors, t-statistic and corresponding (two-sided) p-values.}
   \item{df}{degrees of freedom, a 2-vector, where the first element is the number
      of coefficients and the second element is the number of observations minus
      the number of coefficients.}
   \item{coefCov}{estimated covariance matrix of the coefficients.}
   \item{call*}{the matched call of \code{systemfit}.}
   \item{ols.r.squared*}{OLS \eqn{R^2} value of the entire system.}
   \item{mcelroy.r.squared*}{McElroy's \eqn{R^2} value for the system.}
   \item{iter*}{number of iteration steps (only if the estimation is iterated).}
   \item{control*}{list of control parameters used for the estimation.}
   \item{residCov*}{estimated residual covariance matrix.}
   \item{residCovEst*}{residual covariance matrix used for estimation (only SUR and 3SLS).}
   \item{residCor*}{correlation matrix of the residuals.}
   \item{detResidCov*}{determinant of \code{residCov}.}
   \item{eqnLabel+}{equation label.}
   \item{eqnNo+}{equation number.}
   \item{terms+}{the 'terms' object used for the respective equation.}
   \item{r.squared+}{\eqn{R^2} value of the respective equation.}
   \item{adj.r.squared+}{adjusted \eqn{R^2} value of the respective equation.}
   \item{sigma+}{estimated standard error of the residuals
      of the respective equation.}
   \item{ssr+}{sum of squared residuals of the respective equation.}
   \item{printResidCov*}{argument \code{residCov}.}
   \item{printEquations*}{argument \code{equations}.}
}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com},\cr
   Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de} }

\seealso{\code{\link{systemfit}}, \code{\link{print.systemfit}}}


\examples{
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
inst <- ~ income + farmPrice + trend
system <- list( demand = eqDemand, supply = eqSupply )

## perform OLS on each of the equations in the system
fitols <- systemfit( system, data = Kmenta )

## results of the system
summary( fitols )

## short results of the system
summary( fitols, residCov = FALSE, equations = FALSE )

## results of the first equation
summary( fitols$eq[[1]] )

## results of the second equation
summary( fitols$eq[[2]] )
}

\keyword{models}



