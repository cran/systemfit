
%  $Id: summary.nlsystemfit.system.Rd 154 2006-02-14 14:15:11Z henningsena $

\name{summary.nlsystemfit}
\alias{summary.nlsystemfit.system}
\alias{summary.nlsystemfit.equation}
\title{Summary of nlsystemfit estimation}

\description{
   These functions print a summary of the estimated equation system.
}

\usage{
   \method{summary}{nlsystemfit.system}( object, ... )

   \method{summary}{nlsystemfit.equation}( object, ... )
}

\arguments{
   \item{object}{an object of class \code{nlsystemfit.system} or
      \code{nlsystemfit.equation}.}
   \item{...}{not used by user.}
}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com}}

\seealso{\code{\link{nlsystemfit}}, \code{\link{print.nlsystemfit.system}}}

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

model.3sls <- nlsystemfit( "3SLS", model, start.values, data=ppine,
                                    eqnlabels=labels, inst=inst )
print( model.3sls )
}

\keyword{models}
\keyword{regression}
\keyword{nonlinear}



