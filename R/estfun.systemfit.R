estfun.systemfit <- function ( x, residFit = TRUE, ... ) {
   if( !is.null( x$restrict.matrix ) || !is.null( x$restrict.rhs ) ||
        !is.null( x$restrict.regMat ) ) {
      stop( "returning the estimation function for models with restrictions",
            " has not yet been implemented.")
   }
   
   # residuals
   res <- unlist(  residuals( x ) )
   
   # model matrix
   if( is.null( x$eq[[1]]$inst ) ) {
      mm <- model.matrix( x )
   } else {
      mm <- model.matrix( x, which = "xHat" )
      if( residFit ) {
         res[ !is.na( res ) ] <- res[ !is.na( res ) ] +
            ( model.matrix( x ) - mm ) %*% coef( x )
         #       resid_fit = y - x_fit b
         #       resid = y - x b
         #       resid_fit - resid = - x_fit b + x b
         #       resid_fit = resid + ( x - x_fit ) b
      }
   }
   
   if( sum( !is.na( res ) ) != nrow( mm ) ) {
      stop( "internal error: the number of residuals is not equal to the",
         " number of rows of the model matrix. Please contact the maintainer." )
   }

   if( is.null( x$residCovEst ) ) {
      omegaInvXmat <- mm
   } else {
      omegaInvXmat <- t( .calcXtOmegaInv( xMat = mm, sigma = x$residCovEst, 
         validObsEq = !is.na( residuals( x ) ), invertSigma = TRUE ) )
   }
   
   result <- res[ !is.na( res ) ] * omegaInvXmat
   
   dimnames( result ) <- dimnames( mm )
   
   if( max( abs( colSums( result ) ) ) > 1e-6 ) {
      warning( "the columns of the returned estimating function",
         " do not all sum up to zero,",
         " which indicates that the wrong estimating function is returned" )
   }
   
   return( result )
}