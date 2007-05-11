## this function returns test statistic for
## the hausman test
# The m-statistic is then distributed with k degrees of freedom, where k
# is the rank of the matrix .A generalized inverse is used, as
# recommended by Hausman (1982).

hausman.systemfit <- function( results2sls, results3sls ) {

   result <- list()

   if( is.null( results2sls$bt ) ) {
      result$q <- results2sls$b - results3sls$b
      result$qVar <- results2sls$bcov - results3sls$bcov
   } else {
      result$q <- results2sls$bt - results3sls$bt
      result$qVar <- results2sls$btcov - results3sls$btcov
   }

#    if( min( eigen( hausman$qVar )$values ) < 0 ) {
#       warning( "the matrix V is not 'positive definite'" )
#    }

   result$statistic <- t( result$q ) %*% solve( result$qVar, result$q )
   names( result$statistic ) <- "Hausman"
   result$parameter <- nrow( result$qVar )
   names( result$parameter ) <- "df"
   result$p.value <- 1 - pchisq( result$statistic, result$parameter )
   result$method = paste( "Hausman specification test for consistency of",
      "the 3SLS estimation" )
   result$data.name = results2sls$data.name
   class( result ) <- "htest"
   return( result )
}
