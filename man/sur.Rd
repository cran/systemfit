
% 	$Id: sur.Rd,v 1.3 2002/11/30 02:01:27 hamannj Exp $	

\name{sur}
\alias{sur}
\alias{sur.systemfit}

\title{Seemingly Unrelated Regressions Estimation}

\description{
  Fits a set of structural equations using Seemingly Unrelated
  Regressions. The resulting object is an array of fitting regression
  equations that contain information about the fitting process as well
  as the resulting parameter estimates, standard error estimates and
  covaraince matrix.
}


\usage{
sur.systemfit( eqns, instruments, eqnlabels, data )
}

\arguments{
  \item{eqns}{a list of structural equations to be estimated; a regression
    constant is implied if not explicitly omitted.}
  \item{instruments}{one-sided model formula specifying instrumental
    variables. Instrument variables are not used in this case.}
  \item{data}{an optional data frame containing the variables in the model. 
    By default the variables are taken from the environment from which sur is
    called.}
  \item{eqnlabels}{list of character vectors of names for the equation labels.}
}

\value{
  \code{sur} returns a list of objects of class \code{sur}, with the following components:

  \item{n}{number of observations.}
  \item{p}{number of parameters.}
  \item{coefficients}{parameter estimates.}
  \item{V}{estimated covariance matrix of coefficients.}
  \item{s}{residual standard error.}
  \item{residuals}{vector of residuals.}
  \item{response}{vector of response values.}
  \item{X}{model matrix.}
%  \item{Z}{instrumental-variables matrix.}
  \item{response.name}{name of response variable, or expression evaluating to response.}
  \item{formula}{model formula.}
%  \item{instruments}{one-sided formula for instrumental variables.}

  \item{method}{estimation method for the object, in this case, "SUR".}
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
%  \item{ztzinv}{two-stage instrument regression matrix - ztzinv.}
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
  \item{prediction.limits}{lower and upper 95 percent prediction limits}
  \item{confidence.limits}{lower and upper 95 percent confidence limits}
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

## perform SUR on each of the equations in the system
fitsur <- sur.systemfit( system, inst, labels, kmenta )
print( fitsur )
print( "covariance of residuals" )
print( varcov.systemfit( fitsur ) )

}
\keyword{models}



