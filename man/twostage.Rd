
% 	$Id: twostage.Rd,v 1.2 2002/11/30 02:03:48 hamannj Exp $	

\name{twostage}
\alias{twostage}
\alias{twostage.systemfit}

\title{Two-Stage Least Squares Estimation}

\description{
  Fits a set of structural equations using Two-Stage Least
  Squares. The resulting object is an array of fitting regression
  equations that contain information about the fitting process as well
  as the resulting parameter estimates, standard error estimates and
  covaraince matrix.
}

\usage{
twostage.systemfit( eqns, instruments, eqnlabels, data )
}

\arguments{
  \item{eqns}{a list of structural equations to be estimated; a regression
    constant is implied if not explicitly omitted.}
  \item{instruments}{one-sided model formula specifying instrumental variables.}
  \item{data}{an optional data frame containing the variables in the model. 
    By default the variables are taken from the environment from which twostage is
    called.}
  \item{eqnlabels}{list of character vectors of names for the equation labels.}
}

\value{
  \code{twostage} returns a list of objects of class \code{twostage}, with the following components:

  \item{n}{number of observations.}
  \item{p}{number of parameters.}
  \item{coefficients}{parameter estimates.}
  \item{V}{estimated covariance matrix of coefficients.}
  \item{s}{residual standard error.}
  \item{residuals}{vector of residuals.}
  \item{response}{vector of response values.}
  \item{X}{model matrix.}
  \item{Z}{instrumental-variables matrix.}
  \item{response.name}{name of response variable, or expression evaluating to response.}
  \item{formula}{model formula.}
  \item{instruments}{one-sided formula for instrumental variables.}

  \item{method}{estimation method for the object, in this case, "2SLS".}
  \item{eqnlabel}{the equation label from the labels list for the equation.}			
  \item{formula}{model formula.}
  \item{dfe}{error degrees of freedom.}
  \item{dfm}{model degrees of freedom.}
  \item{model.matrix}{model matrix for the ith equation}
  \item{model.frame}{model frame for the ith equation}
  \item{instruments}{list of instruments for the set of equations.}
  \item{response}{y}
  \item{predicted}{predicted values}
  \item{residuals}{residuals}
  \item{ztzinv}{two-stage instrument regression matrix - ztzinv.}
  \item{v}{v}
  \item{b}{parameter estimates}
  \item{n}{number of observations for the ith equation}
  \item{s}{estimation of sigma}
  \item{sse}{sum of squares}
  \item{mse}{mean squared error} 
  \item{rmse}{square root of mse or root mean squared error}
  \item{se}{estimates standard errors of the parameter estimates.} 
  \item{t}{t values for b}
  \item{p}{p values for b}
  \item{r2}{r2}
  \item{adjr2}{adjusted r-squared}
  \item{covb}{covb}
  \item{se.prediction}{standard error of the predicted values}
  \item{prediction.limts}{lower and upper 95 percent prediction limits}
  \item{confidence.limts}{lower and upper 95 percent confidence limits}
}

\references{    

  Greene, W. H. (1993)
  \emph{Econometric Analysis, Second Edition}, Macmillan.

  Kmenta, J. (1997)
  \emph{Elements of Econometrics, Second Edition}, University of
  Michigan Publishing
  
}

\author{Jeff D. Hamann \email{jeff\_hamann@hamanndonald.com}}

\seealso{\code{\link{ols}},\code{\link{twostage}},\code{\link{sur}}and \code{\link{threestage}}}

\examples{
library( systemfit )

data( kmenta )
attach( kmenta )
demand <- q ~ p + d
supply <- q ~ p + f + a
inst <- ~ d + f + a
labels <- list( "demand", "supply" )
system <- list( demand, supply )

## perform 2SLS on each of the equations in the system
fit2sls <- twostage.systemfit( system, inst, labels, kmenta )
fit3sls <- threestage.systemfit( system, inst, labels, kmenta )

## print the results
print( fit2sls )
print( fit3sls )
print( "covariance of residuals used for estimation (from 2sls)" )
print( varcov.systemfit( fit2sls ) )
print( "covariance of residuals" )
print( varcov.systemfit( fit3sls ) )

## examine the correlation between the predicted values
## of suppy and demand by plotting the correlation over
## the value of q
r12 <- correlation.systemfit( fit3sls, 1, 2 )
plot( q, r12, main="correlation between predictions from supply and demand" )  

## examine the improvement of 3SLS over OLS by computing
## the ratio of the standard errors of the estimates
improve.ratio <- se.ratio.systemfit( fit2sls, fit3sls, 2 )
print( "summary values for the ratio in the std. err. for the predictions" )
print( summary( improve.ratio ) )

}
\keyword{models}



