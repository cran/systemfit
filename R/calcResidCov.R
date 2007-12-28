## Calculate the residual covariance matrix
.calcResidCov <- function( resids, methodResidCov, nObsEq = NULL,
      nCoefEq = NULL, xEq = NULL, diag = FALSE, centered = FALSE,
      useMatrix = FALSE, solvetol = .Machine$double.eps ) {

   eqNames <- NULL
   if( class( resids ) == "data.frame" ) {
      nObsEq <- rep( nrow( resids ), ncol( resids ) )
      eqNames <- names( resids )
      resids <- unlist( resids )
   }
   nEq <- length( nObsEq )
   residi <- list()
   result <- matrix( 0, nEq, nEq )
   for( i in 1:nEq ) {
      residi[[i]] <- resids[ ( 1 + sum(nObsEq[1:i]) - nObsEq[i] ):( sum(nObsEq[1:i]) ) ]
      if( centered ) {
         residi[[i]] <- residi[[i]] - mean( residi[[i]] )
      }
   }
   for( i in 1:nEq ) {
      for( j in ifelse( diag, i, 1 ):ifelse( diag, i, nEq ) ) {
         if( methodResidCov == "noDfCor" ) {
            result[ i, j ] <- sum( residi[[i]] * residi[[j]] ) / nObsEq[i]
         } else if( methodResidCov == "geomean" ) {
            result[ i, j ] <- sum( residi[[i]] * residi[[j]] ) /
               sqrt( ( nObsEq[i] - nCoefEq[i] ) * ( nObsEq[j] - nCoefEq[j] ) )
         } else if( methodResidCov == "Theil" ) {
            #result[ i, j ] <- sum( residi[[i]] * residi[[j]] ) /
            #   ( nObsEq[i] - nCoefEq[i] - nCoefEq[j] + sum( diag(
            #   xEq[[i]] %*% solve( crossprod( xEq[[i]] ), tol=solvetol ) %*%
            #   crossprod( xEq[[i]], xEq[[j]]) %*%
            #   solve( crossprod( xEq[[j]] ), tol=solvetol ) %*%
            #   t( xEq[[j]] ) ) ) )
            result[ i, j ] <- sum( residi[[i]] * residi[[j]] ) /
               ( nObsEq[i] - nCoefEq[i] - nCoefEq[j] + sum( diag(
               solve( crossprod( xEq[[i]] ), tol=solvetol ) %*%
               crossprod( xEq[[i]], xEq[[j]]) %*%
               solve( crossprod( xEq[[j]] ), tol=solvetol ) %*%
               crossprod( xEq[[j]], xEq[[i]] ) ) ) )

         } else if( methodResidCov == "max" ) {
            result[ i, j ] <- sum( residi[[i]] * residi[[j]] ) /
               ( nObsEq[i] - max( nCoefEq[i], nCoefEq[j] ) )
         } else {
            stop( paste( "Argument 'methodResidCov' must be either 'noDfCor',",
                  "'geomean', 'max', or 'Theil'." ) )
         }
      }
   }
   if( !is.null( eqNames ) ) {
      rownames( result ) <- eqNames
      colnames( result ) <- eqNames
   }

   if( useMatrix ){
      result <- as( result, "dspMatrix" )
   }
   return( result )
}

## Calculate Sigma squared
.calcSigma2 <- function( resids, methodResidCov, nObs, nCoef ) {
   if( methodResidCov == "noDfCor" ) {
      result <- sum( resids^2 ) / nObs
   } else if( methodResidCov %in% c( "geomean", "max" ) ){
      result <- sum( resids^2 )/ ( nObs - nCoef )
   } else {
      stop( paste( "Sigma^2 can only be calculated if argument",
         "'methodResidCov' is either 'noDfCor', 'geomean', or 'max'" ) )
   }
}

