waldtest.systemfit <- function( object, R.restr,
   q.restr = rep( 0, nrow( R.restr ) ) ){

   coef <- coef( object )
   vcov <- vcov( object )

   result <- list()

   result$nRestr <- nrow( R.restr )

   result$statistic <- t( R.restr %*% coef - q.restr ) %*%
      solve( R.restr %*% vcov %*% t( R.restr ) ) %*%
      ( R.restr %*% coef - q.restr )

   result$p.value <- 1 - pchisq( result$statistic, result$nRestr )

   class( result ) <- "waldtest.systemfit"
   return( result )
}

print.waldtest.systemfit <- function( x, digits = 4, ... ){
   cat( "\n", "Wald-test for linear parameter restrictions",
      " in equation systems\n", sep = "" )
   cat( "Wald-statistic:", formatC( x$statistic, digits = digits ), "\n" )
   cat( "degrees of freedom:", x$nRestr, "\n" )
   cat( "p-value:", formatC( x$p.value, digits = digits ), "\n\n" )
   invisible( x )
}