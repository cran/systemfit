
library( systemfit )
data( "Kmenta" )

demand <- consump ~ price + income
supply <- consump ~ price + farmPrice + trend
inst <- ~ income + farmPrice + trend
inst1  <- ~ income + farmPrice
instlist <- list( inst1, inst )
labels <- list( "demand", "supply" )
system <- list( demand, supply )
restrm <- matrix(0,1,7)  # restriction matrix "R"
restrm[1,3] <-  1
restrm[1,7] <- -1
restr2m <- matrix(0,2,7)  # restriction matrix "R" 2
restr2q <- matrix(0,2,1)  # restriction vector "q" 2
restr2m[1,3] <-  1
restr2m[1,7] <- -1
restr2m[2,2] <- -1
restr2m[2,5] <-  1
restr2q[2,1] <-  0.5
tc <- matrix(0,7,6)
tc[1,1] <- 1
tc[2,2] <- 1
tc[3,3] <- 1
tc[4,4] <- 1
tc[5,5] <- 1
tc[6,6] <- 1
tc[7,3] <- 1
restr3m <- matrix(0,1,6)  # restriction matrix "R" 2
restr3q <- matrix(0,1,1)  # restriction vector "q" 2
restr3m[1,2] <- -1
restr3m[1,5] <-  1
restr3q[1,1] <-  0.5

# It is not possible to estimate 2SLS with systemfit exactly
# as EViews does, because EViews uses
# rcovformula == 1 for the coefficient covariance matrix and
# rcovformula == 0 for the residual covariance matrix.
# systemfit uses always the same formulas for both calculations.

## *************** 2SLS estimation ************************
## ************ 2SLS estimation (default)*********************
fit2sls1 <- systemfit( "2SLS", system, labels, data = Kmenta, inst = inst )
print( summary( fit2sls1 ) )

## *************** 2SLS estimation (single.eq.sigma=F)*******************
fit2sls1s <- systemfit( "2SLS", system, labels, data = Kmenta, inst = inst,
   single.eq.sigma = FALSE )
print( summary( fit2sls1s ) )

## ********************* 2SLS (probdfsys = TRUE) *****************
fit2sls1p <- systemfit( "2SLS", system, labels, data = Kmenta, inst = inst,
   probdfsys = TRUE )
print( summary( fit2sls1p ) )

## ********************* 2SLS (rcovformula = 0) *****************
fit2sls1r <- systemfit( "2SLS", system, labels, data = Kmenta, inst = inst,
   rcovformula = 0 )
print( summary( fit2sls1r ) )

## *************** 2SLS (rcovformula=0, single.eq.sigma=F) *************
fit2sls1rs <- systemfit( "2SLS", system, labels, data = Kmenta, inst = inst,
   rcovformula = 0, single.eq.sigma = FALSE )
print( summary( fit2sls1rs ) )

## ********************* 2SLS with restriction ********************
## **************** 2SLS with restriction (default)********************
fit2sls2 <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restrm,
   inst = inst )
print( summary( fit2sls2 ) )

## ************* 2SLS with restriction (single.eq.sigma=T) *****************
fit2sls2s <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restrm,
   inst = inst, single.eq.sigma = TRUE )
print( summary( fit2sls2s ) )

## ********************* 2SLS with restriction (probdfsys=T) **************
fit2sls2p <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restrm,
   inst = inst, probdfsys = TRUE )
print( summary( fit2sls2p ) )

## ********************* 2SLS with restriction (rcovformula = 0) **************
fit2sls2r <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restrm,
   inst = inst, rcovformula = 0 )
print( summary( fit2sls2r ) )

## ******** 2SLS with restriction (rcovformula=0, single.eq.sigma=TRUE) *********
fit2sls2rs <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restrm,
   inst = inst, rcovformula = 0, single.eq.sigma = TRUE )
print( summary( fit2sls2rs ) )

## ********************* 2SLS with restriction via TX ******************
## *************** 2SLS with restriction via TX (default )***************
fit2sls3 <- systemfit( "2SLS", system, labels, data = Kmenta, TX = tc,
   inst = inst, rcovformula = 0, probdfsys = TRUE )
print( summary( fit2sls3 ) )

## ********************* 2SLS with restriction via TX (EViews-like) *******
fit2sls3e <- systemfit( "2SLS", system, labels, data = Kmenta, TX = tc,
   inst = inst, rcovformula = 0, probdfsys = TRUE )
print( summary( fit2sls3e ) )

## ***************** 2SLS with 2 restrictions *******************
## ************** 2SLS with 2 restrictions (default) **************
fit2sls4 <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, inst = inst )
print( summary( fit2sls4 ) )

## ************ 2SLS with 2 restrictions (single.eq.sigma=T) **************
fit2sls4s <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, inst = inst, single.eq.sigma = TRUE )
print( summary( fit2sls4s ) )

## ***************** 2SLS with 2 restrictions (probdfsys=T) **************
fit2sls4p <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, inst = inst, probdfsys = TRUE )
print( summary( fit2sls4p ) )

## ***************** 2SLS with 2 restrictions (rcovformula=0) **************
fit2sls4r <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, inst = inst, rcovformula = 0 )
print( summary( fit2sls4r ) )

## ***** 2SLS with 2 restrictions (rcovformula=0, single.eq.sigma=T) *******
fit2sls4rs <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, inst = inst, rcovformula = 0, single.eq.sigma = TRUE )
print( summary( fit2sls4rs ) )

## ************* 2SLS with 2 restrictions via R and TX ******************
## ******** 2SLS with 2 restrictions via R and TX (default) *************
fit2sls5 <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, inst = inst )
print( summary( fit2sls5 ) )

## ******* 2SLS with 2 restrictions via R and TX (single.eq.sigma=T) ******
fit2sls5s <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, inst = inst, single.eq.sigma = TRUE )
print( summary( fit2sls5s ) )

## ********** 2SLS with 2 restrictions via R and TX (probdfsys=T) *******
fit2sls5p <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, inst = inst, probdfsys = TRUE )
print( summary( fit2sls5p ) )

## ************* 2SLS with 2 restrictions via R and TX (rcovformula=0) *********
fit2sls5r <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, inst = inst, rcovformula = 0 )
print( summary( fit2sls5r ) )

## ** 2SLS with 2 restrictions via R and TX (rcovformula=0, single.eq.sigma=T) **
fit2sls5rs <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, inst = inst, rcovformula = 0, single.eq.sigma = TRUE )
print( summary( fit2sls5rs ) )

## *********** 2SLS estimation with different instruments **************
## ******* 2SLS estimation with different instruments (default) *********
fit2slsd1 <- systemfit( "2SLS", system, labels, data = Kmenta, inst = instlist )
print( summary( fit2slsd1 ) )

## *********** 2SLS estimation with different instruments (single.eq.sigma=F)*****
fit2slsd1s <- systemfit( "2SLS", system, labels, data = Kmenta, inst = instlist,
   single.eq.sigma = FALSE )
print( summary( fit2slsd1s ) )

## ********* 2SLS estimation with different instruments (probdfsys=T) *******
fit2slsd1p <- systemfit( "2SLS", system, labels, data = Kmenta, inst = instlist,
   probdfsys = TRUE )
print( summary( fit2slsd1p ) )

## ********* 2SLS estimation with different instruments (rcovformula=0) ******
fit2slsd1r <- systemfit( "2SLS", system, labels, data = Kmenta, inst = instlist,
   rcovformula = 0 )
print( summary( fit2slsd1r ) )

## 2SLS estimation with different instruments (rcovformula=0,single.eq.sigma=F)
fit2slsd1r <- systemfit( "2SLS", system, labels, data = Kmenta, inst = instlist,
   rcovformula = 0, single.eq.sigma = FALSE )
print( summary( fit2slsd1r ) )

## **** 2SLS estimation with different instruments and restriction *******
## ** 2SLS estimation with different instruments and restriction (default) ****
fit2slsd2 <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restrm,
   inst = instlist )
print( summary( fit2slsd2 ) )

## 2SLS estimation with different instruments and restriction (single.eq.sigma=T)
fit2slsd2s <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restrm,
   inst = instlist, single.eq.sigma = TRUE )
print( summary( fit2slsd2s ) )

## **** 2SLS estimation with different instruments and restriction (probdfsys=F)
fit2slsd2p <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restrm,
   inst = instlist, probdfsys = FALSE )
print( summary( fit2slsd2p ) )

## **** 2SLS estimation with different instruments and restriction (rcovformula=0)
fit2slsd2r <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restrm,
   inst = instlist, rcovformula = 0 )
print( summary( fit2slsd2r ) )

## 2SLS estimation with different instr. and restr. (rcovformula=0, single.eq.sigma=T)
fit2slsd2rs <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restrm,
   inst = instlist, rcovformula = 0, single.eq.sigma = TRUE )
print( summary( fit2slsd2rs ) )

## **** 2SLS estimation with different instruments and restriction via TX *
## 2SLS estimation with different instruments and restriction via TX (default)
fit2slsd3 <- systemfit( "2SLS", system, labels, data = Kmenta, TX = tc,
   inst = instlist )
print( summary( fit2slsd3 ) )

## **** 2SLS estimation with different instr. and restr. via TX (rcovformula=0)
fit2slsd3r <- systemfit( "2SLS", system, labels, data = Kmenta, TX = tc,
   inst = instlist, rcovformula = 0 )
print( summary( fit2slsd3r ) )


## *********** estimations with a single regressor ************
fit2slsS1 <- systemfit( "2SLS",
   list( consump ~ price - 1, price ~ consump + trend ),
   data = Kmenta, inst = ~ farmPrice + trend + income )
print( summary( fit2slsS1 ) )
fit2slsS2 <- systemfit( "2SLS",
   list( consump ~ price - 1, consump ~ trend - 1 ),
   data = Kmenta, inst = ~ farmPrice + price + income )
print( summary( fit2slsS2 ) )
fit2slsS3 <- systemfit( "2SLS",
   list( consump ~ trend - 1, price ~ trend - 1 ),
   data = Kmenta, inst = instlist )
print( summary( fit2slsS3 ) )
fit2slsS4 <- systemfit( "2SLS",
   list( consump ~ trend - 1, price ~ trend - 1 ),
   data = Kmenta, inst = ~ farmPrice + trend + income,
   R.restr = matrix( c( 1, -1 ), nrow = 1 ) )
print( summary( fit2slsS4 ) )
fit2slsS5 <- systemfit( "2SLS",
   list( consump ~ 1, price ~ 1 ),
   data = Kmenta, inst = ~ farmPrice )
print( summary( fit2slsS1 ) )


## ****************** residuals **************************
print( residuals( fit2sls1 ) )
print( residuals( fit2sls1$eq[[ 1 ]] ) )

print( residuals( fit2sls2s ) )
print( residuals( fit2sls2s$eq[[ 2 ]] ) )

print( residuals( fit2sls3e ) )
print( residuals( fit2sls3e$eq[[ 1 ]] ) )

print( residuals( fit2sls4r ) )
print( residuals( fit2sls4r$eq[[ 2 ]] ) )

print( residuals( fit2sls5rs ) )
print( residuals( fit2sls5rs$eq[[ 1 ]] ) )

print( residuals( fit2slsd1p ) )
print( residuals( fit2slsd1p$eq[[ 2 ]] ) )

print( residuals( fit2slsd2r ) )
print( residuals( fit2slsd2r$eq[[ 1 ]] ) )


## *********** variance covariance matrix of the coefficients *******
print( round( vcov( fit2sls1s ), digits = 6 ) )
print( round( vcov( fit2sls1s$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fit2sls1r ), digits = 6 ) )
print( round( vcov( fit2sls1r$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fit2sls2p ), digits = 6 ) )
print( round( vcov( fit2sls2p$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fit2sls3 ), digits = 6 ) )
print( round( vcov( fit2sls3$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fit2sls4s ), digits = 6 ) )
print( round( vcov( fit2sls4s$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fit2sls5r ), digits = 6 ) )
print( round( vcov( fit2sls5r$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fit2slsd1 ), digits = 6 ) )
print( round( vcov( fit2slsd1$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fit2slsd2rs ), digits = 6 ) )
print( round( vcov( fit2slsd2rs$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fit2slsd3 ), digits = 6 ) )
print( round( vcov( fit2slsd3$eq[[ 1 ]] ), digits = 6 ) )


## *********** confidence intervals of coefficients *************
print( confint( fit2sls1 ) )
print( confint( fit2sls1$eq[[ 1 ]], level = 0.9 ) )

print( confint( fit2sls2s, level = 0.9 ) )
print( confint( fit2sls2s$eq[[ 2 ]], level = 0.99 ) )

print( confint( fit2sls3e, level = 0.99 ) )
print( confint( fit2sls3e$eq[[ 1 ]], level = 0.5 ) )

print( confint( fit2sls4r, level = 0.5 ) )
print( confint( fit2sls4r$eq[[ 2 ]], level = 0.25 ) )

print( confint( fit2sls5rs, level = 0.25 ) )
print( confint( fit2sls5rs$eq[[ 1 ]], level = 0.975 ) )

print( confint( fit2slsd1p, level = 0.975 ) )
print( confint( fit2slsd1p$eq[[ 2 ]], level = 0.999 ) )

print( confint( fit2slsd2r, level = 0.999 ) )
print( confint( fit2slsd2r$eq[[ 1 ]] ) )


## *********** fitted values *************
print( fitted( fit2sls1, se.fit = TRUE, interval = "prediction" ) )
print( fitted( fit2sls1$eq[[ 1 ]] ) )

print( fitted( fit2sls2s ) )
print( fitted( fit2sls2s$eq[[ 2 ]] ) )

print( fitted( fit2sls3e ) )
print( fitted( fit2sls3e$eq[[ 1 ]] ) )

print( fitted( fit2sls4r ) )
print( fitted( fit2sls4r$eq[[ 2 ]] ) )

print( fitted( fit2sls5rs ) )
print( fitted( fit2sls5rs$eq[[ 1 ]] ) )

print( fitted( fit2slsd1p ) )
print( fitted( fit2slsd1p$eq[[ 2 ]] ) )

print( fitted( fit2slsd2r ) )
print( fitted( fit2slsd2r$eq[[ 1 ]] ) )


## *********** predicted values *************
predictData <- Kmenta
predictData$price <- Kmenta$price * 0.9
predictData$income <- Kmenta$income * 1.1

print( predict( fit2sls1, se.fit = TRUE, interval = "prediction" ) )
print( predict( fit2sls1$eq[[ 1 ]] ) )

print( predict( fit2sls2s, se.pred = TRUE, interval = "confidence",
   level = 0.999, data = predictData ) )
print( predict( fit2sls2s$eq[[ 2 ]] ) )

print( predict( fit2sls3e, se.pred = TRUE, interval = "prediction",
   level = 0.975 ) )
print( predict( fit2sls3e$eq[[ 1 ]] ) )

print( predict( fit2sls4r, se.fit = TRUE, interval = "confidence",
   level = 0.25 ) )
print( predict( fit2sls4r$eq[[ 2 ]] ) )

print( predict( fit2sls5rs, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, data = predictData ) )
print( predict( fit2sls5rs$eq[[ 1 ]] ) )

print( predict( fit2slsd1p, se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.99 ) )
print( predict( fit2slsd1p$eq[[ 2 ]] ) )

print( predict( fit2slsd2r, se.fit = TRUE, interval = "prediction",
   level = 0.9, data = predictData ) )
print( predict( fit2slsd2r$eq[[ 1 ]] ) )

# predict just one observation
smallData <- data.frame( price = 130, income = 150, farmPrice = 120,
   trend = 25, consump = 1 ) ### consump should be removed later!!!!

print( predict( fit2sls1rs, data = smallData ) )
print( predict( fit2sls1rs$eq[[ 1 ]], data = smallData ) )

print( predict( fit2sls2p, se.fit = TRUE, level = 0.9,
   data = smallData ) )
print( predict( fit2sls2p$eq[[ 1 ]], se.pred = TRUE, level = 0.99,
   data = smallData ) )

print( predict( fit2sls3e, interval = "prediction", level = 0.975,
   data = smallData ) )
print( predict( fit2sls3e$eq[[ 1 ]], interval = "confidence", level = 0.8,
   data = smallData ) )

print( predict( fit2sls4r, se.fit = TRUE, interval = "confidence",
   level = 0.999, data = smallData ) )
print( predict( fit2sls4r$eq[[ 2 ]], se.pred = TRUE, interval = "prediction",
   level = 0.75, data = smallData ) )

print( predict( fit2sls5s, se.fit = TRUE, interval = "prediction",
   data = smallData ) )
print( predict( fit2sls5s$eq[[ 1 ]], se.pred = TRUE, interval = "confidence",
   data = smallData ) )

print( predict( fit2slsd3, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, data = smallData ) )
print( predict( fit2slsd3$eq[[ 1 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.25, data = smallData ) )


## ************ correlation of predicted values ***************
print( correlation.systemfit( fit2sls1, 1, 2 ) )

print( correlation.systemfit( fit2sls2s, 2, 1 ) )

print( correlation.systemfit( fit2sls3e, 1, 2 ) )

print( correlation.systemfit( fit2sls4r, 2, 1 ) )

print( correlation.systemfit( fit2sls5rs, 1, 2 ) )

print( correlation.systemfit( fit2slsd1p, 2, 1 ) )

print( correlation.systemfit( fit2slsd2r, 1, 2 ) )


## ************** F tests ****************
# testing first restriction
print( ftest.systemfit( fit2sls1, restrm ) )
print( ftest.systemfit( fit2sls1s, restrm ) )
print( ftest.systemfit( fit2sls1p, restrm ) )
print( ftest.systemfit( fit2sls1r, restrm ) )

# testing second restriction
restrOnly2m <- matrix(0,1,7)
restrOnly2q <- 0.5
restrOnly2m[1,2] <- -1
restrOnly2m[1,5] <-  1
# first restriction not imposed 
print( ftest.systemfit( fit2sls1, restrOnly2m, restrOnly2q ) )
# first restriction imposed
print( ftest.systemfit( fit2sls2, restrOnly2m, restrOnly2q ) )
print( ftest.systemfit( fit2sls2r, restrOnly2m, restrOnly2q ) )
print( ftest.systemfit( fit2sls3, restrOnly2m, restrOnly2q ) )

# testing both of the restrictions
print( ftest.systemfit( fit2sls1, restr2m, restr2q ) )


## ************** Wald tests ****************
# testing first restriction
print( waldtest.systemfit( fit2sls1, restrm ) )
print( waldtest.systemfit( fit2sls1s, restrm ) )
print( waldtest.systemfit( fit2sls1p, restrm ) )
print( waldtest.systemfit( fit2sls1r, restrm ) )

# testing second restriction
# first restriction not imposed
print( waldtest.systemfit( fit2sls1, restrOnly2m, restrOnly2q ) )
# first restriction imposed
print( waldtest.systemfit( fit2sls2, restrOnly2m, restrOnly2q ) )
print( waldtest.systemfit( fit2sls2r, restrOnly2m, restrOnly2q ) )
print( waldtest.systemfit( fit2sls3, restrOnly2m, restrOnly2q ) )

# testing both of the restrictions
print( waldtest.systemfit( fit2sls1, restr2m, restr2q ) )

