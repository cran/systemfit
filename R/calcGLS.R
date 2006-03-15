.calcXtOmegaInv <- function( x, sigma, nObsEq, solvetol = 1e-5, invertSigma = TRUE ){
   nEq <- length( nObsEq )
   if( invertSigma ) {
      sigmaInv <- solve( sigma, tol = solvetol )
   } else {
      sigmaInv <- sigma
   }
   eqSelect <- rep( 0, nrow( x ) )
   for( i in 1:nEq ) {
      eqSelect[ ( sum( nObsEq[ 0:( i - 1 ) ] ) + 1 ):sum( nObsEq[ 1:i ] ) ] <- i
   }
   result <- matrix( 0, nrow = ncol( x ), ncol = nrow( x ) )
   for( i in 1:nEq ) {
      for( j in 1:nEq ) {
         result[ , eqSelect == i ] <- result[ , eqSelect == i ] +
            t( x )[ , eqSelect == j ] * sigmaInv[ i, j ]
      }
   }
   return( result )
}

.calcGLS <- function( x, y = NULL, x2 = x, R.restr = NULL, q.restr = NULL, 
      sigma, nObsEq, solvetol = 1e-5 ){
   xtOmegaInv <- .calcXtOmegaInv( x = x, sigma = sigma, nObsEq = nObsEq, 
      solvetol = solvetol )
   if( is.null( R.restr ) ) {
      result <- solve( xtOmegaInv %*% x2, tol = solvetol )
      if( !is.null( y ) ) {
         result <- result %*% xtOmegaInv %*% y
      }
   } else {
      W <- rbind( cbind( xtOmegaInv %*% x2, t(R.restr) ),
                  cbind( R.restr, matrix(0, nrow(R.restr), nrow(R.restr) )))
      Winv <- solve( W, tol=solvetol )
      if( is.null( y ) ) {
         result <- Winv[ 1:ncol(x), 1:ncol(x) ]
      } else{
         V <- rbind( xtOmegaInv %*% y , q.restr )
         result <- ( Winv %*% V )[1:ncol( x )]     # restricted coefficients
      }
   }
   return( result )
}

