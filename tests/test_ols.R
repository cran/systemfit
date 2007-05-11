
library( systemfit )
data( "Kmenta" )

demand <- consump ~ price + income
supply <- consump ~ price + farmPrice + trend
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

# It is not possible to estimate OLS with systemfit
# exactly as EViews does, because EViews uses
# rcovformula == 1 for the coefficient covariance matrix and
# rcovformula == 0 for the residual covariance matrix, while
# systemfit uses always the same formulas for both calculations.

## *************** OLS estimation ************************
## ********** OLS estimation (default) ********************
fitols1 <- systemfit( "OLS", system, labels, data = Kmenta )
print( summary( fitols1 ) )

## ********** OLS estimation (no single.eq sigma=F) ******************
fitols1s <- systemfit( "OLS", system, labels, data = Kmenta,
   single.eq.sigma = FALSE )
print( summary( fitols1s ) )

## ****************  OLS (probdfsys=T) ***********************
fitols1p <- systemfit( "OLS", system, labels, data = Kmenta,
   probdfsys = TRUE )
print( summary( fitols1p ) )

## ****************  OLS (rcovformula=0) ***********************
fitols1r <- systemfit( "OLS", system, labels, data = Kmenta,
   rcovformula = 0 )
print( summary( fitols1r ) )

## ********  OLS (rcovformula=0, single.eq.sigma=F) ***********
fitols1rs <- systemfit( "OLS", system, labels, data = Kmenta,
   rcovformula = 0, single.eq.sigma = FALSE )
print( summary( fitols1rs ) )

## ****************  OLS (rcovformula=2 ) ***********************
fitols1r <- systemfit( "OLS", system, labels, data = Kmenta,
   rcovformula = 2 )
print( summary( fitols1r ) )

## ****************  OLS (rcovformula=3) ***********************
fitols1r <- systemfit( "OLS", system, labels, data = Kmenta,
   rcovformula = 3 )
print( summary( fitols1r ) )

## ********  OLS (rcovformula=3, single.eq.sigma=F) ***********
fitols1rs <- systemfit( "OLS", system, labels, data = Kmenta,
   rcovformula = 3, single.eq.sigma = FALSE )
print( summary( fitols1rs ) )


## ********* OLS with cross-equation restriction ************
## ****** OLS with cross-equation restriction (default) *********
fitols2 <- systemfit( "OLS", system, labels, data = Kmenta,
   R.restr = restrm )
print( summary( fitols2 ) )

## ****** OLS with cross-equation restriction (single.eq.sigma=T) *******
fitols2s <- systemfit( "OLS", system, labels, data = Kmenta,
   R.restr = restrm, single.eq.sigma = TRUE )
print( summary( fitols2s ) )

## ****** OLS with cross-equation restriction (probdfsys=F) *******
fitols2p <- systemfit( "OLS", system, labels, data = Kmenta,
   R.restr = restrm, probdfsys = FALSE )
print( summary( fitols2p ) )

## ****** OLS with cross-equation restriction (rcovformula=0) *******
fitols2r <- systemfit( "OLS", system, labels, data = Kmenta,
   R.restr = restrm, rcovformula = 0 )
print( summary( fitols2r ) )

## ** OLS with cross-equation restriction (rcovformula=0,single.eq.sigma=T) ***
fitols2rs <- systemfit( "OLS", system, labels, data = Kmenta,
   R.restr = restrm, rcovformula = 0 )
print( summary( fitols2rs ) )

## *** OLS with cross-equation restriction via TX ***
## *** OLS with cross-equation restriction via TX (default) ***
fitols3 <- systemfit( "OLS", system, labels, data = Kmenta, TX = tc )
print( summary( fitols3 ) )

## *** OLS with cross-equation restriction via TX (single.eq.sigma=T) ***
fitols3s <- systemfit( "OLS", system, labels, data = Kmenta,
   TX = tc, single.eq.sigma = TRUE )
print( summary( fitols3s ) )

## *** OLS with cross-equation restriction via TX (probdfsys=F) ***
fitols3p <- systemfit( "OLS", system, labels, data = Kmenta,
   TX = tc, probdfsys = FALSE )
print( summary( fitols3p ) )

## *** OLS with cross-equation restriction via TX (rcovformula=0) ***
fitols3r <- systemfit( "OLS", system, labels, data = Kmenta,
   TX = tc, rcovformula = 0 )
print( summary( fitols3r ) )

## OLS with cross-equation restriction via TX (rcovformula=0,single.eq.sigma=T)
fitols3rs <- systemfit( "OLS", system, labels, data = Kmenta,
   TX = tc, rcovformula = 0, single.eq.sigma = TRUE )
print( summary( fitols3rs ) )

## ********* OLS with 2 cross-equation restrictions ***********
## ********* OLS with 2 cross-equation restrictions (default) ***********
fitols4 <- systemfit( "OLS", system, labels, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q )
print( summary( fitols4 ) )

## ****** OLS with 2 cross-equation restrictions (single.eq.sigma=T) *******
fitols4s <- systemfit( "OLS", system, labels, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, single.eq.sigma = T )
print( summary( fitols4s ) )

## ****** OLS with 2 cross-equation restrictions (probdfsys=F) *******
fitols4p <- systemfit( "OLS", system, labels, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, probdfsys = F )
print( summary( fitols4p ) )

## ****** OLS with 2 cross-equation restrictions (rcovformula=0) *******
fitols4r <- systemfit( "OLS", system, labels, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, rcovformula = 0 )
print( summary( fitols4r ) )

## OLS with 2 cross-equation restrictions (rcovformula=0, single.eq.sigma=T) *
fitols4rs <- systemfit( "OLS", system, labels, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, rcovformula = 0, single.eq.sigma = T )
print( summary( fitols4rs ) )

## ***** OLS with 2 cross-equation restrictions via R and TX ****
## ***** OLS with 2 cross-equation restrictions via R and TX (default) ****
fitols5 <- systemfit( "OLS", system, labels, data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, rcovformula = 0)
print( summary( fitols5 ) )

## ***** OLS with 2 cross-equation restrictions via R and TX (single.eq.sigma=T) ****
fitols5s <- systemfit( "OLS", system, labels, data = Kmenta,R.restr = restr3m,
   q.restr = restr3q, TX = tc, single.eq.sigma = T )
print( summary( fitols5s ) )

## ***** OLS with 2 cross-equation restrictions via R and TX (probdfsys=F) ****
fitols5p <- systemfit( "OLS", system, labels, data = Kmenta,R.restr = restr3m,
   q.restr = restr3q, TX = tc, probdfsys = F )
print( summary( fitols5p ) )

## ***** OLS with 2 cross-equation restrictions via R and TX (rcovformula=0) ****
fitols5r <- systemfit( "OLS", system, labels, data = Kmenta,R.restr = restr3m,
   q.restr = restr3q, TX = tc, rcovformula = 0 )
print( summary( fitols5r ) )

## OLS with 2 cross-equation restr. via R and TX (rcovformula=0,single.eq.sigma=T)
fitols5rs <- systemfit( "OLS", system, labels, data = Kmenta,R.restr = restr3m,
   q.restr = restr3q, TX = tc, rcovformula = 0, single.eq.sigma = T )
print( summary( fitols5rs ) )


## ****************** residuals **************************
print( residuals( fitols1p ) )
print( residuals( fitols1p$eq[[ 2 ]] ) )

print( residuals( fitols2r ) )
print( residuals( fitols2r$eq[[ 1 ]] ) )

print( residuals( fitols3s ) )
print( residuals( fitols3s$eq[[ 2 ]] ) )

print( residuals( fitols4rs ) )
print( residuals( fitols4rs$eq[[ 1 ]] ) )

print( residuals( fitols5 ) )
print( residuals( fitols5$eq[[ 2 ]] ) )


## *********** variance covariance matrix of the coefficients *******
print( round( vcov( fitols1rs ), digits = 6 ) )
print( round( vcov( fitols1rs$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitols2s ), digits = 6 ) )
print( round( vcov( fitols2s$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitols3p ), digits = 6 ) )
print( round( vcov( fitols3p$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitols4r ), digits = 6 ) )
print( round( vcov( fitols4r$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitols5 ), digits = 6 ) )
print( round( vcov( fitols5$eq[[ 2 ]] ), digits = 6 ) )


## *********** confidence intervals of coefficients *************
print( confint( fitols1p ) )
print( confint( fitols1p$eq[[ 2 ]], level = 0.9 ) )

print( confint( fitols2r, level = 0.9 ) )
print( confint( fitols2r$eq[[ 1 ]], level = 0.99 ) )

print( confint( fitols3s, level = 0.99 ) )
print( confint( fitols3s$eq[[ 2 ]], level = 0.5 ) )

print( confint( fitols4rs, level = 0.5 ) )
print( confint( fitols4rs$eq[[ 1 ]], level = 0.25 ) )

print( confint( fitols5, level = 0.25 ) )
print( confint( fitols5$eq[[ 2 ]], level = 0.999 ) )

print( confint( fitols3p, level = 0.999 ) )
print( confint( fitols3p$eq[[ 1 ]] ) )


## *********** fitted values *************
print( fitted( fitols1p ) )
print( fitted( fitols1p$eq[[ 2 ]] ) )

print( fitted( fitols2r ) )
print( fitted( fitols2r$eq[[ 1 ]] ) )

print( fitted( fitols3s ) )
print( fitted( fitols3s$eq[[ 2 ]] ) )

print( fitted( fitols4rs ) )
print( fitted( fitols4rs$eq[[ 1 ]] ) )

print( fitted( fitols5 ) )
print( fitted( fitols5$eq[[ 2 ]] ) )


## *********** predicted values *************
predictData <- Kmenta
predictData$price <- Kmenta$price * 0.9
predictData$income <- Kmenta$income * 1.1

print( predict( fitols1p, se.fit = TRUE, interval = "prediction" ) )
print( predict( fitols1p$eq[[ 2 ]] ) )

print( predict( fitols2r, se.pred = TRUE, interval = "confidence",
   level = 0.999, data = predictData ) )
print( predict( fitols2r$eq[[ 1 ]] ) )

print( predict( fitols3s, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, data = predictData ) )
print( predict( fitols3s$eq[[ 2 ]] ) )

print( predict( fitols4rs, se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.99 ) )
print( predict( fitols4rs$eq[[ 1 ]] ) )

print( predict( fitols5, se.fit = TRUE, interval = "prediction",
   level = 0.9, data = predictData ) )
print( predict( fitols5$eq[[ 2 ]] ) )

# predict just one observation
smallData <- data.frame( price = 130, income = 150, farmPrice = 120,
   trend = 25, consump = 1 ) ### consump should be removed later!!!!

print( predict( fitols1p, data = smallData ) )
print( predict( fitols1p$eq[[ 1 ]], data = smallData ) )

print( predict( fitols2r, se.fit = TRUE, level = 0.9,
   data = smallData ) )
print( predict( fitols2r$eq[[ 1 ]], se.pred = TRUE, level = 0.99,
   data = smallData ) )

print( predict( fitols3s, interval = "prediction", level = 0.975,
   data = smallData ) )
print( predict( fitols3s$eq[[ 1 ]], interval = "confidence", level = 0.8,
   data = smallData ) )

print( predict( fitols4rs, se.fit = TRUE, interval = "confidence",
   level = 0.999, data = smallData ) )
print( predict( fitols4rs$eq[[ 2 ]], se.pred = TRUE, interval = "prediction",
   level = 0.75, data = smallData ) )

print( predict( fitols5, se.fit = TRUE, interval = "prediction",
   data = smallData ) )
print( predict( fitols5$eq[[ 1 ]], se.pred = TRUE, interval = "confidence",
   data = smallData ) )

print( predict( fitols5rs, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, data = smallData ) )
print( predict( fitols5rs$eq[[ 1 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.25, data = smallData ) )


## ************ correlation of predicted values ***************
print( correlation.systemfit( fitols1p, 1, 2 ) )

print( correlation.systemfit( fitols2r, 2, 1 ) )

print( correlation.systemfit( fitols3s, 1, 2 ) )

print( correlation.systemfit( fitols4rs, 2, 1 ) )

print( correlation.systemfit( fitols5, 1, 2 ) )


## ************** F tests ****************
# testing first restriction
print( ftest.systemfit( fitols1, restrm ) )
print( ftest.systemfit( fitols1s, restrm ) )
print( ftest.systemfit( fitols1p, restrm ) )
print( ftest.systemfit( fitols1r, restrm ) )

# testing second restriction
restrOnly2m <- matrix(0,1,7)
restrOnly2q <- 0.5
restrOnly2m[1,2] <- -1
restrOnly2m[1,5] <-  1
# first restriction not imposed 
print( ftest.systemfit( fitols1, restrOnly2m, restrOnly2q ) )
# first restriction imposed
print( ftest.systemfit( fitols2, restrOnly2m, restrOnly2q ) )
print( ftest.systemfit( fitols3, restrOnly2m, restrOnly2q ) )

# testing both of the restrictions
print( ftest.systemfit( fitols1, restr2m, restr2q ) )


## ************** Wald tests ****************
# testing first restriction
print( waldtest.systemfit( fitols1, restrm ) )
print( waldtest.systemfit( fitols1s, restrm ) )
print( waldtest.systemfit( fitols1p, restrm ) )
print( waldtest.systemfit( fitols1r, restrm ) )

# testing second restriction
# first restriction not imposed
print( waldtest.systemfit( fitols1, restrOnly2m, restrOnly2q ) )
# first restriction imposed
print( waldtest.systemfit( fitols2, restrOnly2m, restrOnly2q ) )
print( waldtest.systemfit( fitols3, restrOnly2m, restrOnly2q ) )

# testing both of the restrictions
print( waldtest.systemfit( fitols1, restr2m, restr2q ) )
