.systemfitPanel <- function( formula, data, pooled ) {

   if( class( data )[1] != "pdata.frame" ) {
      stop( "argument 'data' must be of class 'pdata.frame'" )
   }
   eqnVar <- attributes( data )$indexes$id
   timeVar <- attributes( data )$indexes$time

   result <- list()

   data[[ eqnVar ]] <- gsub( " |_", ".", data[[ eqnVar ]] )
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
         newVar <- paste( eqn, var, sep = "_" )
         newVar <- make.names( newVar )
         wideData[[ newVar ]] <- NA
         wideData[ make.names( eqnData[ , timeVar ] ), newVar ] <-
            eqnData[ , var ]
         endogVar <- gsub( var, newVar, endogVar )
         exogVar <- gsub( var, newVar, exogVar )
      }
      eqnSystem[[ eqnNo ]] <- as.formula( paste( endogVar, "~", exogVar ) )
   }
   names( eqnSystem ) <- eqnLabels
   #reshape( data, idvar=timevar, timevar=eqnVar,direction="wide")

   restrict.regMat <- NULL
   if( pooled ) {
      for( eqnNo in 1:nEqn ) {
         restrict.regMat <- rbind( restrict.regMat, diag( 1, nRegressors ) )
      }
   }

   result$eqnSystem <- eqnSystem
   result$wideData  <- wideData
   result$restrict.regMat <- restrict.regMat

   return( result )
}