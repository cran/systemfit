library( "systemfit" )
library( "sandwich" )
options( warn = 1 )

data( "KleinI" )
eqConsump  <- consump ~ corpProf + corpProfLag + wages
eqInvest   <- invest ~ corpProf + corpProfLag + capitalLag
eqPrivWage <- privWage ~ gnp + gnpLag + trend
inst <- ~ govExp + taxes + govWage + trend + capitalLag + corpProfLag + gnpLag
system <- list( Consumption = eqConsump, Investment = eqInvest,
   PrivateWages = eqPrivWage )
restrict <- c( "Consumption_corpProf + Investment_capitalLag = 0" )
restrict2 <- c( restrict, "Consumption_corpProfLag - PrivateWages_trend = 0" )

for( dataNo in 1:5 ) {
   # set some values of some variables to NA
   if( dataNo == 2 ) {
      KleinI$gnpLag[ 7 ] <- NA
   } else if( dataNo == 3 ) {
      KleinI$wages[ 10 ] <- NA
   } else if( dataNo == 4 ) {
      KleinI$corpProf[ 13 ] <- NA
   } else if( dataNo == 5 ) {
      KleinI$invest[ 16 ] <- NA
   }

for( methodNo in 1:5 ) {
   method <- c( "OLS", "2SLS", "SUR", "3SLS", "3SLS" )[ methodNo ]
   maxit <- ifelse( methodNo == 5, 500, 1 )

   cat( "> \n> # ", ifelse( maxit == 1, "", "I" ), method, "\n", sep = "" )
   if( method %in% c( "OLS", "WLS", "SUR" ) ) {
      kleinModel <- systemfit( system, method = method, data = KleinI,
         methodResidCov = ifelse( method == "OLS", "geomean", "noDfCor" ),
         maxit = maxit )
   } else {
      kleinModel <- systemfit( system, method = method, data = KleinI,
         inst = inst, methodResidCov = "noDfCor", maxit = maxit )
   }
   cat( "> summary\n" )
   print( summary( kleinModel ) )
   cat( "> residuals\n" )
   print( residuals( kleinModel ) )
   cat( "> fitted\n" )
   print( fitted( kleinModel ) )
   cat( "> predict\n" )
   print( predict( kleinModel, se.fit = TRUE,
      interval = ifelse( methodNo %in% c( 1, 4 ), "prediction",  "confidence" ),
      useDfSys = methodNo %in% c( 1, 3, 5 ) ) )
   cat( "> model.frame\n" )
   if( methodNo == 1 ) {
      mfOls <- model.frame( kleinModel )
      print( mfOls )
   } else if( methodNo == 2 ) {
      mf2sls <- model.frame( kleinModel )
      print( mf2sls )
   cat( "> Frames of instrumental variables\n" )
      for( i in 1:3 ){
         print( kleinModel$eq[[ i ]]$modelInst )
      }
   } else if( methodNo == 3 ) {
      print( all.equal( mfOls, model.frame( kleinModel ) ) )
   } else {
      print( all.equal( mf2sls, model.frame( kleinModel ) ) )
   }
   cat( "> model.matrix\n" )
   if( methodNo == 1 ) {
      mmOls <- model.matrix( kleinModel )
      print( mmOls )
   } else {
      print( all.equal( mmOls, model.matrix( kleinModel ) ) )
   }
   if( methodNo == 2 ) {
      cat( "> matrix of instrumental variables\n" )
      print( model.matrix( kleinModel, which = "z" ) )
      cat( "> matrix of fitted regressors\n" )
      print( model.matrix( kleinModel, which = "xHat" ) )
   }
   cat( "> nobs\n" )
   print( nobs( kleinModel ) )
   cat( "> linearHypothesis\n" )
   print( linearHypothesis( kleinModel, restrict ) )
   print( linearHypothesis( kleinModel, restrict, test = "F" ) )
   print( linearHypothesis( kleinModel, restrict, test = "Chisq" ) )
   print( linearHypothesis( kleinModel, restrict2 ) )
   print( linearHypothesis( kleinModel, restrict2, test = "F" ) )
   print( linearHypothesis( kleinModel, restrict2, test = "Chisq" ) )
   cat( "> logLik\n" )
   print( logLik( kleinModel ) )
   
   cat( "Estimating function\n" )
   print( estfun( kleinModel ) )
   print( colSums( estfun( kleinModel ) ) )

   cat( "> Bread\n" )
   print( bread( kleinModel ) )
}
}