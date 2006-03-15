systemfitClassic <- function( method, formula, eqnVar, timeVar, data,
   pooled = FALSE, ... ) {

   eqnLabels <- levels( as.factor( data[[ eqnVar ]] ) )
   nEqn <- length( eqnLabels )
   timeLabels <- levels( as.factor( data[[ timeVar ]] ) )
   nRegressors <- ncol( model.matrix( formula, data ) )
   eqnSystem <- list()

   wideData <- data.frame( time = timeLabels )
   rownames( wideData ) <- make.names( timeLabels )
   for( eqnNo in 1:nEqn ) {
      eqn <- eqnLabels[ eqnNo ]
      endogVar <- formula[2]
      exogVar <- formula[3]
      eqnData <- data[ data[[ eqnVar ]] == eqn, ]
      for( var in all.vars( formula ) ) {
         newVar <- paste( var, eqn, sep = "." )
         newVar <- make.names( newVar )
         wideData[[ newVar ]] <- NA
         wideData[ make.names( eqnData[ , timeVar ] ), newVar ] <-
            eqnData[ , var ]
         endogVar <- gsub( var, newVar, endogVar )
         exogVar <- gsub( var, newVar, exogVar )
      }
      eqnSystem[[ eqnNo ]] <- as.formula( paste( endogVar, "~", exogVar ) )
   }
   #reshape( data, idvar=timevar, timevar=eqnVar,direction="wide")

   TX <- NULL
   if( pooled ) {
      for( eqnNo in 1:nEqn ) {
         TX <- rbind( TX, diag( 1, nRegressors ) )
      }
   }

   result <- systemfit( method = method, eqns = eqnSystem,
      eqnlabels = eqnLabels, data = wideData, TX = TX, ... )

   return( result )
}