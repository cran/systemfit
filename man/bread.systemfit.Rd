\name{bread.systemfit}
\alias{bread.systemfit}

\title{Bread for Sandwiches}

\description{
   Extract the estimator for the bread of sandwhiches
   (see \code{\link[sandwich]{bread}}).
}

\usage{
\method{bread}{systemfit}( obj, ... )
}

\arguments{
   \item{obj}{an object of class \code{systemfit}.}
   \item{\dots}{further arguments (currently ignored).}
}

\value{
   Quadratic symmetric matrix,
   which is an estimator for the expectation of the negative derivative 
   of the estimating function (see \code{\link{estfun.systemfit}}).
}

\section{Warnings}{
   The \pkg{sandwich} package must be loaded before this method
   can be used.
   
   This method might not be suitable for specific formulas 
   for 3SLS estimations in case of unbalanced systems
   or different instruments for different equations.
}

\author{
   Arne Henningsen
}

\seealso{\code{\link[sandwich]{bread}}, \code{\link{systemfit}}.}

\examples{
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
system <- list( demand = eqDemand, supply = eqSupply )
inst <- ~ income + farmPrice + trend

## OLS estimation
fitols <- systemfit( system, "OLS", data = Kmenta )

## obtain the bread
library( "sandwich" )
bread( fitols )

## this is only true for OLS models
all.equal( bread( fitols ),
   solve( crossprod( model.matrix( fitols ) ) / 40 ) )


## 2SLS estimation
fit2sls <- systemfit( system, "2SLS", inst = inst, data = Kmenta )

## obtain the bread
bread( fit2sls )

## this is only true for 2SLS models
all.equal( bread( fit2sls ),
   solve( crossprod( model.matrix( fit2sls, which = "xHat" ) ) / 40 ) )


## iterated SUR estimation
fitsur <- systemfit( system, "SUR", data = Kmenta, maxit = 100 )

## obtain the bread
bread( fitsur )

## this should be true for SUR and WLS models
all.equal( bread( fitsur ),
   solve( t( model.matrix( fitsur ) ) \%*\%
      (  ( solve( fitsur$residCovEst ) \%x\% diag( nrow( Kmenta ) ) ) \%*\% 
      model.matrix( fitsur ) ) / 40 ), check.attributes = FALSE )


## 3SLS estimation
fit3sls <- systemfit( system, "3SLS", inst = inst, data = Kmenta )

## obtain the bread
bread( fit3sls )

## this should be true for 3SLS and W2SLS models
all.equal( bread( fit3sls ),
   solve( t( model.matrix( fit3sls, which = "xHat" ) ) \%*\%
      (  ( solve( fit3sls$residCovEst ) \%x\% diag( nrow( Kmenta ) ) ) \%*\% 
      model.matrix( fit3sls, which = "xHat" ) ) / 40 ), check.attributes = FALSE )
}

\keyword{methods}
