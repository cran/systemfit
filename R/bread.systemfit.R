bread.systemfit <- function ( x, ... ) {
   if( !is.null( x$restrict.matrix ) || !is.null( x$restrict.rhs ) ||
      !is.null( x$restrict.regMat ) ) {
      stop( "returning the 'bread' for models with restrictions",
            " has not yet been implemented.")
   }
   
   # model matrix
   if( is.null( x$eq[[1]]$inst ) ) {
      mm <- model.matrix( x )
   } else {
      mm <- model.matrix( x, which = "xHat" )
   }
   
   if( is.null( x$residCovEst ) ) {
      omegaInvXmat <- mm
   } else {
      omegaInvXmat <- t( .calcXtOmegaInv( xMat = mm, sigma = x$residCovEst, 
         validObsEq = !is.na( residuals( x ) ), invertSigma = TRUE ) )
   }
   
   result <- solve( crossprod( mm, omegaInvXmat ) / nrow( mm ) )
   
   return( result )
}