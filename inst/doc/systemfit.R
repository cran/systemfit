## ----include=FALSE------------------------------------------------------------
library(knitr)
opts_chunk$set(
engine='R'
)

## ----echo=FALSE---------------------------------------------------------------
options( prompt = "R> ", ctinue = "+  " )

## ----eval=FALSE---------------------------------------------------------------
#  sigmaInv <- solve( residCov )
#  xtOmegaInv <- crossprod( xMat, kronecker( sigmaInv, Diagonal( nObs ) ) )
#  coef <- solve( xtOmegaInv %*% xMat, xtOmegaInv %*%  yVec )

## ----message=FALSE------------------------------------------------------------
library( "systemfit" )
data( "Kmenta" )
attach( Kmenta )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
eqSystem <- list( demand = eqDemand, supply = eqSupply )
fitols <- systemfit( eqSystem )
print( fitols )

## -----------------------------------------------------------------------------
fitsur <- systemfit( eqSystem, method = "SUR" )

## -----------------------------------------------------------------------------
fit3sls  <- systemfit( eqSystem, method = "3SLS",
   inst = ~ income + farmPrice + trend )
fit3sls2 <- systemfit( eqSystem, method = "3SLS",
   inst = list( ~ farmPrice + trend, ~ income + farmPrice + trend ) )

## -----------------------------------------------------------------------------
fitsur <- systemfit( eqSystem, method = "SUR", data = Kmenta )

## ----tidy=FALSE---------------------------------------------------------------
restrict <- "demand_price + supply_farmPrice = 0"
fitsurRmat <- systemfit(eqSystem, method = "SUR",
   restrict.matrix = restrict)

## ----tidy=FALSE---------------------------------------------------------------
Rmat <- matrix(0, nrow = 1, ncol = 7)
Rmat[1, 2] <- 1
Rmat[1, 6] <- 1
qvec <- c(0)
fitsurRmatNum <- systemfit(eqSystem, method = "SUR",
   restrict.matrix = Rmat, restrict.rhs = qvec)

## ----tidy=FALSE---------------------------------------------------------------
modRegMat <- matrix(0, nrow = 7, ncol = 6)
modRegMat[1:5, 1:5] <- diag(5)
modRegMat[6, 2] <- -1
modRegMat[7, 6] <- 1
fitsurRegMat <- systemfit(eqSystem, method = "SUR",
   restrict.regMat = modRegMat)

## -----------------------------------------------------------------------------
all.equal( coef( fitsurRmat ), coef( fitsurRmatNum ) )
all.equal( coef( fitsurRmat ), coef( fitsurRegMat ) )

## -----------------------------------------------------------------------------
fitsurit <- systemfit( eqSystem, method = "SUR", maxiter = 500 )

## -----------------------------------------------------------------------------
summary( fitsur )

## -----------------------------------------------------------------------------
summary( fitsur, residCov = FALSE, equations = FALSE )

## ----message=FALSE,eval = requireNamespace("plm", quietly = TRUE)-------------
### this code chunk is evaluated only if the 'plm' package is available
data( "GrunfeldGreene" )
library( "plm" )
GGPanel <- pdata.frame( GrunfeldGreene, c( "firm", "year" ) )
greeneSur <- systemfit( invest ~ value + capital, method = "SUR",
   data = GGPanel )

## ----eval = requireNamespace("plm", quietly = TRUE)---------------------------
### this code chunk is evaluated only if the 'plm' package is available
greeneSurPooled <- systemfit( invest ~ value + capital, method = "SUR",
   data = GGPanel, pooled = TRUE )

## -----------------------------------------------------------------------------
linearHypothesis( fitsur, Rmat, qvec, test = "FT" )

linearHypothesis( fitsur, Rmat, qvec, test = "F" )

linearHypothesis( fitsur, Rmat, qvec, test = "Chisq" )

lrtest( fitsurRmat, fitsur )

## -----------------------------------------------------------------------------
fit2sls  <- systemfit( eqSystem, method = "2SLS",
   inst = ~ income + farmPrice + trend, data = Kmenta )
fit3sls <- systemfit( eqSystem, method = "3SLS",
   inst = ~ income + farmPrice + trend, data = Kmenta )
hausman.systemfit( fit2sls, fit3sls )

## ----echo=FALSE---------------------------------------------------------------
hausmantest <- hausman.systemfit( fit2sls, fit3sls )

## ----echo=FALSE---------------------------------------------------------------
library( "systemfit" )

## ----results='hide'-----------------------------------------------------------
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
inst <- ~ income + farmPrice + trend
system <- list( demand = eqDemand, supply = eqSupply )

## ----results='hide'-----------------------------------------------------------
fitOls <- systemfit( system, data = Kmenta )
round( coef( summary( fitOls ) ), digits = 4 )

## ----results='hide'-----------------------------------------------------------
fit2sls <- systemfit( system, method = "2SLS", inst = inst, data = Kmenta )
round( coef( summary( fit2sls ) ), digits = 4 )

## ----results='hide'-----------------------------------------------------------
fit3sls <- systemfit( system, method = "3SLS", inst = inst, data = Kmenta )
round( coef( summary( fit3sls ) ), digits = 4 )

## ----results='hide'-----------------------------------------------------------
fitI3sls <- systemfit( system, method = "3SLS", inst = inst, data = Kmenta,
   maxit = 250 )
round( coef( summary( fitI3sls ) ), digits = 4 )

## -----------------------------------------------------------------------------
data( "KleinI" )
eqConsump  <- consump ~ corpProf + corpProfLag + wages
eqInvest   <- invest ~ corpProf + corpProfLag + capitalLag
eqPrivWage <- privWage ~ gnp + gnpLag + trend
inst <- ~ govExp + taxes + govWage + trend + capitalLag + corpProfLag + gnpLag
system <- list( Consumption = eqConsump, Investment = eqInvest,
   PrivateWages = eqPrivWage )

## ----results='hide'-----------------------------------------------------------
kleinOls <- systemfit( system, data = KleinI )
round( coef( summary( kleinOls ) ), digits = 3 )

## ----results='hide'-----------------------------------------------------------
klein2sls <- systemfit( system, method = "2SLS", inst = inst, data = KleinI,
   methodResidCov = "noDfCor" )
round( coef( summary( klein2sls ) ), digits = 3 )

## ----results='hide'-----------------------------------------------------------
klein3sls <- systemfit( system, method = "3SLS", inst = inst, data = KleinI,
   methodResidCov = "noDfCor" )
round( coef( summary( klein3sls ) ), digits = 3 )

## ----results='hide'-----------------------------------------------------------
kleinI3sls <- systemfit( system, method = "3SLS", inst = inst, data = KleinI,
   methodResidCov = "noDfCor", maxit = 500 )
round( coef( summary( kleinI3sls ) ), digits = 3 )

## ----message=FALSE,eval = requireNamespace("plm", quietly = TRUE)-------------
### this code chunk is evaluated only if the 'plm' package is available
data( "GrunfeldGreene" )
library( "plm" )
GGPanel <- pdata.frame( GrunfeldGreene, c( "firm", "year" ) )
formulaGrunfeld <- invest ~ value + capital

## ----results='hide',eval = requireNamespace("plm", quietly = TRUE)------------
### this code chunk is evaluated only if the 'plm' package is available
greeneOls <- systemfit( formulaGrunfeld,
   data = GGPanel )
round( coef( summary( greeneOls ) ), digits = 4 )
round( sapply( greeneOls$eq, function(x){return(summary(x)$ssr/20)} ), digits = 3 )

## ----results='hide',eval = requireNamespace("plm", quietly = TRUE)------------
### this code chunk is evaluated only if the 'plm' package is available
greeneOlsPooled <- systemfit( formulaGrunfeld,
   data = GGPanel, pooled = TRUE )
round( coef( summary( greeneOlsPooled$eq[[1]] ) ), digits = 4 ) #$
sum( sapply( greeneOlsPooled$eq, function(x){return(summary(x)$ssr)}) )/100

## ----results='hide',eval = requireNamespace("plm", quietly = TRUE)------------
### this code chunk is evaluated only if the 'plm' package is available
greeneSur <- systemfit( formulaGrunfeld, method = "SUR",
   data = GGPanel, methodResidCov = "noDfCor" )
round( coef( summary( greeneSur ) ), digits = 4 )
round( greeneSur$residCov, digits = 3 ) #$
round( summary( greeneSur )$residCor, digits = 3 ) #$

## ----results='hide',eval = requireNamespace("plm", quietly = TRUE)------------
### this code chunk is evaluated only if the 'plm' package is available
greeneSurPooled <- systemfit( formulaGrunfeld, method = "SUR",
   data = GGPanel, pooled = TRUE, methodResidCov = "noDfCor",
   residCovWeighted = TRUE )
round( coef( summary( greeneSurPooled$eq[[1]] ) ), digits = 4 ) #$

round( greeneSurPooled$residCov, digits = 3 ) #$
round( cov( residuals( greeneSurPooled ) ), digits = 3 )
round( summary( greeneSurPooled )$residCor, digits = 3 ) #$

## ----eval = requireNamespace("plm", quietly = TRUE)---------------------------
### this code chunk is evaluated only if the 'plm' package is available
GrunfeldTheil <- subset( GrunfeldGreene,
   firm %in% c( "General Electric", "Westinghouse" ) )
GTPanel <- pdata.frame( GrunfeldTheil, c( "firm", "year" ) )
formulaGrunfeld <- invest ~ value + capital

## ----results='hide',eval = requireNamespace("plm", quietly = TRUE)------------
### this code chunk is evaluated only if the 'plm' package is available
theilOls <- systemfit( formulaGrunfeld,
   data = GTPanel )
round( coef( summary( theilOls ) ), digits = 3 )

## ----results='hide',eval = requireNamespace("plm", quietly = TRUE)------------
### this code chunk is evaluated only if the 'plm' package is available
theilSur <- systemfit( formulaGrunfeld, method = "SUR",
   data = GTPanel, methodResidCov = "noDfCor" )
round( coef( summary( theilSur ) ), digits = 3 )

## ----results='hide',eval = requireNamespace("plm", quietly = TRUE)------------
### this code chunk is evaluated only if the 'plm' package is available
RMatrix <- rbind( c( 0, 1, 0, 0, -1, 0 ), c( 0, 0, 1, 0, 0, -1 ) )
linearHypothesis( theilSur, RMatrix )

## ----results='hide',tidy=FALSE,eval = requireNamespace("plm", quietly = TRUE)----
### this code chunk is evaluated only if the 'plm' package is available
theilSurRestr <- systemfit(formulaGrunfeld, method = "SUR",
   data = GTPanel, methodResidCov = "noDfCor", restrict.matrix = RMatrix,
   residCovRestricted = FALSE)
round(coef(summary(theilSurRestr)), digits = 3)

## ----message=FALSE------------------------------------------------------------
library( "systemfit" )
set.seed( 1 )
systemfitModel <- createSystemfitModel( nEq = 8, nReg = 10, nObs = 750 )
system.time(
   fitMatrix <- systemfit( systemfitModel$formula, method = "SUR",
      data = systemfitModel$data )
)
system.time(
   fitTrad <- systemfit( systemfitModel$formula, method = "SUR",
      data = systemfitModel$data, useMatrix = FALSE )
)
all.equal( fitMatrix, fitTrad )

## -----------------------------------------------------------------------------
smallModel <- createSystemfitModel( nEq = 3, nReg = 4, nObs = 50 )
system.time(
   fitSmallMatrix <- systemfit( smallModel$formula, method = "SUR",
      data = smallModel$data, maxit = 500 )
)
system.time(
   fitSmallTrad <- systemfit( smallModel$formula, method = "SUR",
      data = smallModel$data, maxit = 500, useMatrix = FALSE )
)
all.equal( fitSmallMatrix, fitSmallTrad )

## ----echo=FALSE----------------------------------------------------------
options( width = 75 )

## ----message=FALSE,eval = requireNamespace("sem", quietly = TRUE)--------
### this code chunk is evaluated only if the 'sem' package is available
library( "sem" )
library( "systemfit" )
data( "KleinI" )

## ----tidy=FALSE,eval = requireNamespace("sem", quietly = TRUE)-----------
### this code chunk is evaluated only if the 'sem' package is available
limlRam <- matrix(c(
   "consump  <-  corpProf",    "consump_corpProf",    NA,
   "consump  <-  corpProfLag", "consump_corpProfLag", NA,
   "consump  <-  wages",       "consump_wages",       NA,
   "invest   <-  corpProf",    "invest_corpProf",     NA,
   "invest   <-  corpProfLag", "invest_corpProfLag",  NA,
   "invest   <-  capitalLag",  "invest_capitalLag",   NA,
   "privWage <-  gnp",         "privWage_gnp",        NA,
   "privWage <-  gnpLag",      "privWage_gnpLag",     NA,
   "privWage <-  trend",       "privWage_trend",      NA,
   "consump  <-> consump",     "s11", NA,
   "privWage <-> privWage",    "s22", NA,
   "invest   <-> invest",      "s33", NA),
   ncol = 3, byrow = TRUE)
class(limlRam) <- "mod"
exogVar <- c("corpProf", "corpProfLag", "wages", "capitalLag", "trend",
   "gnp", "gnpLag")
endogVar <- c("consump", "invest", "privWage")
allVar <- c(exogVar, endogVar)

limlResult <- sem(model = limlRam, S = cov(KleinI[ -1, allVar ]),
   N = (nrow(KleinI) - 1), fixed.x = exogVar)
print(limlResult)

## ----tidy=FALSE----------------------------------------------------------
eqConsump <- consump ~ corpProf + corpProfLag + wages
eqInvest <- invest ~ corpProf + corpProfLag + capitalLag
eqPrivWage <- privWage ~ gnp + gnpLag + trend
system <- list(consump = eqConsump, invest = eqInvest,
   privWage = eqPrivWage)
olsResult <- systemfit(system, data = KleinI)
print(olsResult)

## ----tidy=FALSE,eval = requireNamespace("sem", quietly = TRUE)-----------
### this code chunk is evaluated only if the 'sem' package is available
fimlRam <- rbind(limlRam,
   c("consump  <-> invest",   "s12", NA),
   c("consump  <-> privWage", "s13", NA),
   c("privWage <-> invest",   "s23", NA))
class(fimlRam) <- "mod"

fimlResult <- sem(model = fimlRam, S = cov(KleinI[ -1, allVar ]),
   N = (nrow(KleinI) - 1), fixed.x = exogVar)
print(fimlResult)

## ------------------------------------------------------------------------
surResult <- systemfit( system, method = "SUR", data = KleinI,
   methodResidCov = "noDfCor", maxit = 500 )
print( surResult )

