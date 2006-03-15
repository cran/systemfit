## Likelihood Ratio Test
lrtest.systemfit <- function( resultc, resultu ) {
  lrtest <- list()
  if( resultc$method %in% c( "SUR", "WSUR" ) &
      resultu$method %in% c( "SUR", "WSUR" ) ) {
    n   <- resultu$eq[[1]]$n
    lrtest$nRestr  <- resultu$ki - resultc$ki
    if(resultc$rcovformula != resultu$rcovformula) {
      stop( paste( "both estimations must use the same formula to calculate",
                   "the residual covariance matrix!" ) )
    }
    if(resultc$rcovformula == 0) {
      lrtest$statistic  <- n * ( log( resultc$drcov ) - log( resultu$drcov ) )
    } else {
      residc <- array(resultc$resids,c(n,resultc$g))
      residu <- array(resultu$resids,c(n,resultu$g))
      lrtest$statistic <- n * ( log( det( (t(residc) %*% residc)) ) -
                         log( det( (t(residu) %*% residu))))
    }
    lrtest$p.value <- 1 - pchisq( lrtest$statistic, lrtest$nRestr )
  }
  class( lrtest ) <- "lrtest.systemfit"
  return( lrtest )
}

print.lrtest.systemfit <- function( x, digits = 4, ... ){
   cat( "\n", "Likelihood-Ratio-test for parameter restrictions",
      " in equation systems\n", sep = "" )
   cat( "LR-statistic:", formatC( x$statistic, digits = digits ), "\n" )
   cat( "degrees of freedom:", x$nRestr, "\n" )
   cat( "p-value:", formatC( x$p.value, digits = digits ), "\n\n" )
   invisible( x )
}

