logLik.systemfit <- function( object, ... ){

   resid <- residuals( object )
   residCov <- .calcResidCov( resid, "noDfCor" )
   nObsPerEq <- nrow( resid )
   nEq <- ncol( resid )

   result <- - ( nObsPerEq / 2 ) * ( nEq * ( 1 + log( 2 * pi ) ) +
      log( det( residCov ) ) )

   nCoef <- nEq * nObsPerEq - df.residual( object )
   if( object$method %in% c( "OLS", "2SLS" ) ){
      nSigma <- 1
   } else if( object$method %in% c( "WLS", "W2SLS" ) ){
      nSigma <- nEq
   } else if( object$method %in% c( "SUR", "3SLS" ) ){
      nSigma <- nEq * ( nEq + 1 ) / 2
   } else {
      stop( "internal error: unknown estimation method '", object$method, "'" )
   }

   attributes( result )$nobs <- nEq * nObsPerEq
   attributes( result )$df <- nCoef + nSigma
   class( result ) <- "logLik"

   return( result )
}