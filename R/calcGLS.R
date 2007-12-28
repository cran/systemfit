.calcXtOmegaInv <- function( xMat, sigma, nObsEq, invertSigma = TRUE,
      useMatrix = FALSE, warnMatrix = TRUE, solvetol = 1e-5 ){

   nEq <- length( nObsEq )

   if( useMatrix && warnMatrix ){
      if( class( sigma ) != "dspMatrix" ){
         warning( "class of 'sigma' is '", class( sigma ),
            "', but it should be 'dspMatrix'" )
      }
      if( class( xMat ) != "dgCMatrix" ){
         warning( "class of 'xMat' is '", class( xMat ),
            "', but it should be 'dgCMatrix'" )
      }
   }

   if( invertSigma ) {
      sigmaInv <- solve( sigma, tol = solvetol )
   } else {
      sigmaInv <- sigma
   }

   if( useMatrix ){
      result <- crossprod( xMat, suppressWarnings(
         kronecker( sigmaInv, Diagonal( nObsEq[ 1 ] ) ) ) )
   } else {
      eqSelect <- rep( 0, nrow( xMat ) )
      for( i in 1:nEq ) {
         eqSelect[ ( sum( nObsEq[ 0:( i - 1 ) ] ) + 1 ):sum( nObsEq[ 1:i ] ) ] <- i
      }
      result <- matrix( 0, nrow = ncol( xMat ), ncol = nrow( xMat ) )
      for( i in 1:nEq ) {
         for( j in 1:nEq ) {
            result[ , eqSelect == i ] <- result[ , eqSelect == i ] +
               t( xMat )[ , eqSelect == j ] * sigmaInv[ i, j ]
         }
      }
   }

   return( result )
}

.calcGLS <- function( xMat, yVec = NULL, xMat2 = xMat, R.restr = NULL,
      q.restr = NULL, sigma, nObsEq, useMatrix = TRUE, warnMatrix = TRUE,
      solvetol = 1e-5 ){

   if( useMatrix && warnMatrix ){
      if( class( xMat ) != "dgCMatrix" ){
         warning( "class of 'xMat' is '", class( xMat ),
            "', but it should be 'dgCMatrix'" )
      }
      if( class( xMat2 ) != "dgCMatrix" ){
         warning( "class of 'xMat2' is '", class( xMat2 ),
            "', but it should be 'dgCMatrix'" )
      }
      if( class( sigma ) != "dspMatrix" ){
         warning( "class of 'sigma' is '", class( sigma ),
            "', but it should be 'dspMatrix'" )
      }
   }

   xtOmegaInv <- .calcXtOmegaInv( xMat = xMat, sigma = sigma, nObsEq = nObsEq,
      useMatrix = useMatrix, solvetol = solvetol )
   if( is.null( R.restr ) ) {
      if( is.null( yVec ) ) {
         result <- solve( xtOmegaInv %*% xMat2, tol = solvetol )
      } else {
         result <- as.numeric( solve( xtOmegaInv %*% xMat2, xtOmegaInv %*% yVec,
            tol = solvetol ) )
      }
   } else {
      W <- rbind2( cbind2( xtOmegaInv %*% xMat2, t(R.restr) ),
                  cbind2( R.restr, matrix(0, nrow(R.restr), nrow(R.restr) )))
      if( is.null( yVec ) ) {
         result <- as.matrix(
            solve( W, tol=solvetol )[ 1:ncol(xMat), 1:ncol(xMat) ] )
      } else{
         V <- c( as.numeric( xtOmegaInv %*% yVec ), q.restr )
         result <- solve( W, V, tol=solvetol )[ 1:ncol( xMat ) ]
      }
   }
   return( result )
}

