## Calculate the residual covariance matrix
.calcRCov <- function( resids, rcovformula, nObsEq = NULL, nCoefEq = NULL, xEq = NULL,
      diag = FALSE, centered = FALSE, solvetol = .Machine$double.eps ) {

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
         if( rcovformula == 0 ) {
            result[ i, j ] <- sum( residi[[i]] * residi[[j]] ) / nObsEq[i]
         } else if( rcovformula == 1 || rcovformula == "geomean" ) {
            result[ i, j ] <- sum( residi[[i]] * residi[[j]] ) /
               sqrt( ( nObsEq[i] - nCoefEq[i] ) * ( nObsEq[j] - nCoefEq[j] ) )
         } else if( rcovformula == 2 || rcovformula == "Theil" ) {
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

         } else if( rcovformula == 3 || rcovformula == "max" ) {
            result[ i, j ] <- sum( residi[[i]] * residi[[j]] ) /
               ( nObsEq[i] - max( nCoefEq[i], nCoefEq[j] ) )
         } else {
            stop( paste( "Argument 'rcovformula' must be either 0, 1,",
                  "'geomean', 2, 'Theil', 3 or 'max'." ) )
         }
      }
   }
   if( !is.null( eqNames ) ) {
      rownames( result ) <- eqNames
      colnames( result ) <- eqNames
   }
   return( result )
}

## Calculate Sigma squared
.calcSigma2 <- function( resids, rcovformula, nObs, nCoef ) {
   if( rcovformula == 0 ) {
      result <- sum( resids^2 ) / nObs
   } else if( rcovformula == 1 || rcovformula == "geomean" ||
      rcovformula == 3 || rcovformula == "max") {
      result <- sum( resids^2 )/ ( nObs - nCoef )
   } else {
      stop( paste( "Sigma^2 can only be calculated if argument",
         "'rcovformula' is either 0, 1, 'geomean', 3 or 'max'" ) )
   }
}

