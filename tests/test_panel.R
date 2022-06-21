library( systemfit )
if(requireNamespace( 'plm', quietly = TRUE ) ) {
   library( plm )
   options( digits = 3 )
   useMatrix <- FALSE
}

## Repeating the OLS and SUR estimations in Theil (1971, pp. 295, 300)
if(requireNamespace( 'plm', quietly = TRUE ) ) {
   data( "GrunfeldGreene" )
   GrunfeldTheil <- subset( GrunfeldGreene,
      firm %in% c( "General Electric", "Westinghouse" ) )
   GrunfeldTheil <- pdata.frame( GrunfeldTheil, c( "firm", "year" ) )
   formulaGrunfeld <- invest ~ value + capital
}

# OLS
if(requireNamespace( 'plm', quietly = TRUE ) ) {
   theilOls <- systemfit( formulaGrunfeld, "OLS",
      data = GrunfeldTheil, useMatrix = useMatrix )
   print( theilOls )
   print( summary( theilOls ) )
   print( summary( theilOls, useDfSys = TRUE, residCov = FALSE,
      equations = FALSE ) )
   print( summary( theilOls, equations = FALSE ) )
   print( coef( theilOls ) )
   print( coef( summary(theilOls ) ) )
   print( vcov( theilOls ) )
   print( residuals( theilOls ) )
   print( confint( theilOls ) )
   print( fitted(theilOls ) )
   print( logLik( theilOls ) )
   print( logLik( theilOls, residCovDiag = TRUE ) )
   print( nobs( theilOls ) )
   print( model.frame( theilOls ) )
   print( model.matrix( theilOls ) )
   print( formula( theilOls ) )
   print( formula( theilOls$eq[[ 1 ]] ) )
   print( terms( theilOls ) )
   print( terms( theilOls$eq[[ 1 ]] ) )
}

# SUR
if(requireNamespace( 'plm', quietly = TRUE ) ) {
   theilSur <- systemfit( formulaGrunfeld, "SUR",
      data = GrunfeldTheil, methodResidCov = "noDfCor", useMatrix = useMatrix )
   print( theilSur )
   print( summary( theilSur ) )
   print( summary( theilSur, useDfSys = TRUE, equations = FALSE ) )
   print( summary( theilSur, residCov = FALSE, equations = FALSE ) )
   print( coef( theilSur ) )
   print( coef( summary( theilSur ) ) )
   print( vcov( theilSur ) )
   print( residuals( theilSur ) )
   print( confint( theilSur ) )
   print( fitted( theilSur ) )
   print( logLik( theilSur ) )
   print( logLik( theilSur, residCovDiag = TRUE ) )
   print( nobs( theilSur ) )
   print( model.frame( theilSur ) )
   print( model.matrix( theilSur ) )
   print( formula( theilSur ) )
   print( formula( theilSur$eq[[ 2 ]] ) )
   print( terms( theilSur ) )
   print( terms( theilSur$eq[[ 2 ]] ) )
}

## Repeating the OLS and SUR estimations in Greene (2003, pp. 351)
if(requireNamespace( 'plm', quietly = TRUE ) ) {
   GrunfeldGreene <- pdata.frame( GrunfeldGreene, c( "firm", "year" ) )
   formulaGrunfeld <- invest ~ value + capital
}

# OLS
if(requireNamespace( 'plm', quietly = TRUE ) ) {
   greeneOls <- systemfit( formulaGrunfeld, "OLS",
      data = GrunfeldGreene, useMatrix = useMatrix )
   print( greeneOls )
   print( summary( greeneOls ) )
   print( summary( greeneOls, useDfSys = TRUE, equations = FALSE ) )
   print( summary( greeneOls, residCov = FALSE ) )
   print( sapply( greeneOls$eq, function(x){return(summary(x)$ssr/20)} ) ) # sigma^2
   print( coef( greeneOls ) )
   print( coef( summary( greeneOls ) ) )
   print( vcov( greeneOls ) )
   print( residuals( greeneOls ) )
   print( confint(greeneOls  ) )
   print( fitted( greeneOls ) )
   print( logLik( greeneOls ) )
   print( logLik( greeneOls, residCovDiag = TRUE ) )
   print( nobs( greeneOls ) )
   print( model.frame( greeneOls ) )
   print( model.matrix( greeneOls ) )
   print( formula( greeneOls ) )
   print( formula( greeneOls$eq[[ 2 ]] ) )
   print( terms( greeneOls ) )
   print( terms( greeneOls$eq[[ 2 ]] ) )
}

# OLS Pooled
if(requireNamespace( 'plm', quietly = TRUE ) ) {
   greeneOlsPooled <- systemfit( formulaGrunfeld, "OLS",
      data = GrunfeldGreene, pooled = TRUE, useMatrix = useMatrix )
   print( greeneOlsPooled )
   print( summary( greeneOlsPooled ) )
   print( summary( greeneOlsPooled, useDfSys = FALSE, residCov = FALSE ) )
   print( summary( greeneOlsPooled, residCov = FALSE, equations = FALSE ) )
   print( sum( sapply( greeneOlsPooled$eq, function(x){return(summary(x)$ssr)}) )/97 ) # sigma^2
   print( coef( greeneOlsPooled ) )
   print( coef( greeneOlsPooled, modified.regMat = TRUE ) )
   print( coef( summary( greeneOlsPooled ) ) )
   print( coef( summary( greeneOlsPooled ), modified.regMat = TRUE ) )
   print( vcov( greeneOlsPooled ) )
   print( vcov( greeneOlsPooled, modified.regMat = TRUE ) )
   print( residuals( greeneOlsPooled ) )
   print( confint( greeneOlsPooled ) )
   print( fitted( greeneOlsPooled ) )
   print( logLik( greeneOlsPooled ) )
   print( logLik( greeneOlsPooled, residCovDiag = TRUE ) )
   print( nobs( greeneOlsPooled ) )
   print( model.frame( greeneOlsPooled ) )
   print( model.matrix( greeneOlsPooled ) )
   print( formula( greeneOlsPooled ) )
   print( formula( greeneOlsPooled$eq[[ 1 ]] ) )
   print( terms( greeneOlsPooled ) )
   print( terms( greeneOlsPooled$eq[[ 1 ]] ) )
}

# SUR
if(requireNamespace( 'plm', quietly = TRUE ) ) {
   greeneSur <- systemfit( formulaGrunfeld, "SUR",
      data = GrunfeldGreene, methodResidCov = "noDfCor", useMatrix = useMatrix )
   print( greeneSur )
   print( summary( greeneSur ) )
   print( summary( greeneSur, useDfSys = TRUE, residCov = FALSE ) )
   print( summary( greeneSur, equations = FALSE ) )
   print( coef( greeneSur ) )
   print( coef( summary( greeneSur ) ) )
   print( vcov( greeneSur ) )
   print( residuals( greeneSur ) )
   print( confint( greeneSur ) )
   print( fitted( greeneSur ) )
   print( logLik( greeneSur ) )
   print( logLik( greeneSur, residCovDiag = TRUE ) )
   print( nobs( greeneSur ) )
   print( model.frame( greeneSur ) )
   print( model.matrix( greeneSur ) )
   print( formula( greeneSur ) )
   print( formula( greeneSur$eq[[ 1 ]] ) )
   print( terms( greeneSur ) )
   print( terms( greeneSur$eq[[ 1 ]] ) )
}

# SUR Pooled
if(requireNamespace( 'plm', quietly = TRUE ) ) {
   greeneSurPooled <- systemfit( formulaGrunfeld, "SUR",
      data = GrunfeldGreene, pooled = TRUE, methodResidCov = "noDfCor",
      residCovWeighted = TRUE, useMatrix = useMatrix )
   print( greeneSurPooled )
   print( summary( greeneSurPooled ) )
   print( summary( greeneSurPooled, useDfSys = FALSE, equations = FALSE ) )
   print( summary( greeneSurPooled, residCov = FALSE, equations = FALSE ) )
   print( coef( greeneSurPooled ) )
   print( coef( greeneSurPooled, modified.regMat = TRUE ) )
   print( coef( summary( greeneSurPooled ) ) )
   print( coef( summary( greeneSurPooled ), modified.regMat = TRUE ) )
   print( vcov( greeneSurPooled ) )
   print( vcov( greeneSurPooled, modified.regMat = TRUE ) )
   print( residuals( greeneSurPooled ) )
   print( confint( greeneSurPooled ) )
   print( fitted( greeneSurPooled ) )
   print( logLik( greeneSurPooled ) )
   print( logLik( greeneSurPooled, residCovDiag = TRUE ) )
   print( nobs( greeneSurPooled ) )
   print( model.frame( greeneSurPooled ) )
   print( model.matrix( greeneSurPooled ) )
   print( formula( greeneSurPooled ) )
   print( formula( greeneSurPooled$eq[[ 1 ]] ) )
   print( terms( greeneSurPooled ) )
   print( terms( greeneSurPooled$eq[[ 1 ]] ) )
}


#########  IV estimation  #######################
###  2SLS  ###
# instruments = explanatory variables  ->  2SLS estimates = OLS estimates
if(requireNamespace( 'plm', quietly = TRUE ) ) {
   greene2sls <- systemfit( formulaGrunfeld, inst = ~ value + capital, "2SLS",
      data = GrunfeldGreene, useMatrix = useMatrix )
   print( greene2sls )
   print( summary( greene2sls ) )
   print( all.equal( coef( summary( greene2sls ) ), coef( summary( greeneOls ) ) ) )
   print( all.equal( greene2sls[ -c(1,2,6) ], greeneOls[ -c(1,2,6) ] ) )
   for( i in 1:length( greene2sls$eq ) ) {
      print( all.equal( greene2sls$eq[[i]][ -c(3,15:17) ], 
         greeneOls$eq[[i]][-3] ) )
   }
}
# 'real' IV/2SLS estimation
if(requireNamespace( 'plm', quietly = TRUE ) ) {
   greene2slsR <- systemfit( invest ~ capital, inst = ~ value, "2SLS",
      data = GrunfeldGreene, useMatrix = useMatrix )
   print( greene2slsR )
   print( summary( greene2slsR ) )
}

###  2SLS, pooled  ###
# instruments = explanatory variables  ->  2SLS estimates = OLS estimates
if(requireNamespace( 'plm', quietly = TRUE ) ) {
   greene2slsPooled <- systemfit( formulaGrunfeld, inst = ~ value + capital, "2SLS",
      data = GrunfeldGreene, pooled = TRUE, useMatrix = useMatrix )
   print( greene2slsPooled )
   print( summary( greene2slsPooled ) )
   print( all.equal( coef( summary( greene2slsPooled ) ),
      coef( summary( greeneOlsPooled ) ) ) )
   print( all.equal( greene2slsPooled[ -c(1,2,6) ], greeneOlsPooled[ -c(1,2,6) ] ) )
   for( i in 1:length( greene2slsPooled$eq ) ) {
      print( all.equal( greene2slsPooled$eq[[i]][ -c(3,15:17) ], 
         greeneOlsPooled$eq[[i]][-3] ) )
   }
}
# 'real' IV/2SLS estimation
if(requireNamespace( 'plm', quietly = TRUE ) ) {
   greene2slsRPooled <- systemfit( invest ~ capital, inst = ~ value, "2SLS",
      data = GrunfeldGreene, pooled = TRUE, useMatrix = useMatrix )
   print( greene2slsRPooled )
   print( summary( greene2slsRPooled ) )
}

###  3SLS  ###
# instruments = explanatory variables  ->  3SLS estimates = SUR estimates
if(requireNamespace( 'plm', quietly = TRUE ) ) {
   greene3sls <- systemfit( formulaGrunfeld, inst = ~ value + capital, "3SLS",
      data = GrunfeldGreene, useMatrix = useMatrix, methodResidCov = "noDfCor" )
   print( greene3sls )
   print( summary( greene3sls ) )
   print( all.equal( coef( summary( greene3sls ) ), coef( summary( greeneSur ) ) ) )
   print( all.equal( greene3sls[ -c(1,2,7) ], greeneSur[ -c(1,2,7) ] ) )
   for( i in 1:length( greene3sls$eq ) ) {
      print( all.equal( greene3sls$eq[[i]][ -c(3,15:17) ], 
         greeneSur$eq[[i]][-3] ) )
   }
}
# 'real' IV/3SLS estimation
if(requireNamespace( 'plm', quietly = TRUE ) ) {
   greene3slsR <- systemfit( invest ~ capital, inst = ~ value, "3SLS",
      data = GrunfeldGreene, useMatrix = useMatrix )
   print( greene3slsR )
   print( summary( greene3slsR ) )
}

###  3SLS, Pooled  ###
# instruments = explanatory variables  ->  3SLS estimates = SUR estimates
if(requireNamespace( 'plm', quietly = TRUE ) ) {
   greene3slsPooled <- systemfit( formulaGrunfeld, inst = ~ capital + value, "3SLS",
      data = GrunfeldGreene, pooled = TRUE, useMatrix = useMatrix, 
      residCovWeighted = TRUE, methodResidCov = "noDfCor" )
   print( greene3slsPooled )
   print( summary( greene3slsPooled ) )
   print( all.equal( coef( summary( greene3slsPooled ) ),
      coef( summary( greeneSurPooled ) ) ) )
   print( all.equal( greene3slsPooled[ -c(1,2,7) ], greeneSurPooled[ -c(1,2,7) ] ) )
   for( i in 1:length( greene3slsPooled$eq ) ) {
      print( all.equal( greene3slsPooled$eq[[i]][ -c(3,15:17) ], 
         greeneSurPooled$eq[[i]][-3] ) )
   }
}
# 'real' IV/3SLS estimation
if(requireNamespace( 'plm', quietly = TRUE ) ) {
   greene3slsRPooled <- systemfit( invest ~ capital, inst = ~ value, "3SLS",
      data = GrunfeldGreene, useMatrix = useMatrix )
   print( greene3slsRPooled )
   print( summary( greene3slsRPooled ) )
}


## **************** estfun ************************
library( "sandwich" )

if(requireNamespace( 'plm', quietly = TRUE ) ) {
   print( estfun( theilOls ) )
   print( round( colSums( estfun( theilOls ) ), digits = 7 ) )
   
   print( estfun( theilSur ) )
   print( round( colSums( estfun( theilSur ) ), digits = 7 ) )
   
   print( estfun( greeneOls ) )
   print( round( colSums( estfun( greeneOls ) ), digits = 7 ) )
   
   print( try( estfun( greeneOlsPooled ) ) )
   
   print( estfun( greeneSur ) )
   print( round( colSums( estfun( greeneSur ) ), digits = 7 ) )
   
   print( try( estfun( greeneSurPooled ) ) )
}

## **************** bread ************************
if(requireNamespace( 'plm', quietly = TRUE ) ) {
   print( bread( theilOls ) )
   
   print( bread( theilSur ) )
   
   print( bread( greeneOls ) )
   
   print( try( bread( greeneOlsPooled ) ) )
   
   print( bread( greeneSur ) )
   
   print( try( bread( greeneSurPooled ) ) )
}
