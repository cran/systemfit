
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


## ********************* W2SLS *****************
fitw2sls1 <- systemfit( "W2SLS", system, labels, data = Kmenta, inst = inst )
print( summary( fitw2sls1 ) )

## ********************* W2SLS (EViews-like) *****************
fitw2sls1e <- systemfit( "W2SLS", system, labels, data = Kmenta, inst = inst,
   rcovformula = 0, probdfsys = TRUE )
print( summary( fitw2sls1e ) )

## ********************* W2SLS with restriction *******************
fitw2sls2 <- systemfit( "W2SLS", system, labels, data = Kmenta, R.restr = restrm,
   inst = inst )
print( summary( fitw2sls2 ) )

## ********************* W2SLS with restriction (EViews-like) **************
fitw2sls2e <- systemfit( "W2SLS", system, labels, data = Kmenta, R.restr = restrm,
   inst = inst, rcovformula = 0, probdfsys = TRUE )
print( summary( fitw2sls2e ) )

## ********************* W2SLS with restriction via TX *******************
fitw2sls3 <- systemfit( "W2SLS", system, labels, data = Kmenta, TX = tc, inst = inst )
print( summary( fitw2sls3 ) )

## ********************* W2SLS with restriction via TX (EViews-like) **************
fitw2sls3e <- systemfit( "W2SLS", system, labels, data = Kmenta, TX = tc, inst = inst,
   rcovformula = 0, probdfsys = TRUE )
print( summary( fitw2sls3e ) )

## ***************** W2SLS with 2 restrictions ********************
fitw2sls4 <- systemfit( "W2SLS", system, labels, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, inst = inst )
print( summary( fitw2sls4 ) )

## ***************** W2SLS with 2 restrictions (EViews-like) **************
fitw2sls4e <- systemfit( "W2SLS", system, labels, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, inst = inst, rcovformula = 0, probdfsys = TRUE )
print( summary( fitw2sls4e ) )

## ***************** W2SLS with 2 restrictions via R and TX ******************
fitw2sls5 <- systemfit( "W2SLS", system, labels, data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, inst = inst )
print( summary( fitw2sls5 ) )

## ***************** W2SLS with 2 restrictions via R and TX (EViews-like) **************
fitw2sls5e <- systemfit( "W2SLS", system, labels, data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, inst = inst, rcovformula = 0, probdfsys = TRUE )
print( summary( fitw2sls5e ) )

## ****** 2SLS estimation with different instruments **********************
fitw2slsd1 <- systemfit( "W2SLS", system, labels, data = Kmenta, inst = instlist )
print( summary( fitw2slsd1 ) )

## ****** 2SLS estimation with different instruments (EViews-like)******************
fitw2slsd1e <- systemfit( "W2SLS", system, labels, data = Kmenta, inst = instlist,
   rcovformula = 0, probdfsys = TRUE )
print( summary( fitw2slsd1e ) )

## **** W2SLS estimation with different instruments and restriction ********
fitw2slsd2 <- systemfit( "W2SLS", system, labels, data = Kmenta, R.restr = restrm,
   inst = instlist )
print( summary( fitw2slsd2 ) )

## **** W2SLS estimation with different instruments and restriction (EViews-like)*
fitw2slsd2e <- systemfit( "W2SLS", system, labels, data = Kmenta, R.restr = restrm,
   inst = instlist, rcovformula = 0, probdfsys = TRUE )
print( summary( fitw2slsd2e ) )

## ** W2SLS estimation with different instruments and restriction via TX ****
fitw2slsd3 <- systemfit( "W2SLS", system, labels, data = Kmenta, TX = tc,
   inst = instlist)
print( summary( fitw2slsd3 ) )

## W2SLS estimation with different instruments and restriction via TX (EViews-like)
fitw2slsd3e <- systemfit( "W2SLS", system, labels, data = Kmenta, TX = tc,
   inst = instlist, rcovformula = 0, probdfsys = TRUE )
print( summary( fitw2slsd3e ) )


## *********** estimations with a single regressor ************
fitw2slsS1 <- systemfit( "W2SLS",
   list( consump ~ price - 1, price ~ consump + trend ),
   data = Kmenta, inst = ~ farmPrice + trend + income )
print( summary( fitw2slsS1 ) )
fitw2slsS2 <- systemfit( "W2SLS",
   list( consump ~ price - 1, consump ~ trend - 1 ),
   data = Kmenta, inst = ~ farmPrice + price + income )
print( summary( fitw2slsS2 ) )
fitw2slsS3 <- systemfit( "W2SLS",
   list( consump ~ trend - 1, price ~ trend - 1 ),
   data = Kmenta, inst = instlist )
print( summary( fitw2slsS3 ) )
fitw2slsS4 <- systemfit( "W2SLS",
   list( consump ~ trend - 1, price ~ trend - 1 ),
   data = Kmenta, inst = ~ farmPrice + trend + income,
   R.restr = matrix( c( 1, -1 ), nrow = 1 ) )
print( summary( fitw2slsS4 ) )
fitw2slsS5 <- systemfit( "W2SLS",
   list( consump ~ 1, price ~ 1 ),
   data = Kmenta, inst = instlist )
print( summary( fitw2slsS5 ) )


## ****************** residuals **************************
print( residuals( fitw2sls1e ) )
print( residuals( fitw2sls1e$eq[[ 1 ]] ) )

print( residuals( fitw2sls2 ) )
print( residuals( fitw2sls2$eq[[ 2 ]] ) )

print( residuals( fitw2sls3 ) )
print( residuals( fitw2sls3$eq[[ 1 ]] ) )

print( residuals( fitw2sls4e ) )
print( residuals( fitw2sls4e$eq[[ 2 ]] ) )

print( residuals( fitw2sls5 ) )
print( residuals( fitw2sls5$eq[[ 1 ]] ) )

print( residuals( fitw2slsd1 ) )
print( residuals( fitw2slsd1$eq[[ 2 ]] ) )

print( residuals( fitw2slsd2e ) )
print( residuals( fitw2slsd2e$eq[[ 1 ]] ) )

print( residuals( fitw2slsd3e ) )
print( residuals( fitw2slsd3e$eq[[ 2 ]] ) )


## *********** variance covariance matrix of the coefficients *******
print( round( vcov( fitw2sls1e ), digits = 6 ) )
print( round( vcov( fitw2sls1e$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitw2sls2 ), digits = 6 ) )
print( round( vcov( fitw2sls2$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitw2sls3e ), digits = 6 ) )
print( round( vcov( fitw2sls3e$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitw2sls4 ), digits = 6 ) )
print( round( vcov( fitw2sls4$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitw2sls5 ), digits = 6 ) )
print( round( vcov( fitw2sls5$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitw2slsd1 ), digits = 6 ) )
print( round( vcov( fitw2slsd1$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitw2slsd2e ), digits = 6 ) )
print( round( vcov( fitw2slsd2e$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitw2slsd3 ), digits = 6 ) )
print( round( vcov( fitw2slsd3$eq[[ 1 ]] ), digits = 6 ) )


## *********** confidence intervals of coefficients *************
print( confint( fitw2sls1e ) )
print( confint( fitw2sls1e$eq[[ 1 ]], level = 0.9 ) )

print( confint( fitw2sls2, level = 0.9 ) )
print( confint( fitw2sls2$eq[[ 2 ]], level = 0.99 ) )

print( confint( fitw2sls3, level = 0.99 ) )
print( confint( fitw2sls3$eq[[ 1 ]], level = 0.5 ) )

print( confint( fitw2sls4e, level = 0.5 ) )
print( confint( fitw2sls4e$eq[[ 2 ]], level = 0.25 ) )

print( confint( fitw2sls5, level = 0.25 ) )
print( confint( fitw2sls5$eq[[ 1 ]], level = 0.975 ) )

print( confint( fitw2slsd1, level = 0.975 ) )
print( confint( fitw2slsd1$eq[[ 2 ]], level = 0.999 ) )

print( confint( fitw2slsd2e, level = 0.999 ) )
print( confint( fitw2slsd2e$eq[[ 1 ]], level = 0.01 ) )

print( confint( fitw2slsd3e, level = 0.01 ) )
print( confint( fitw2slsd3e$eq[[ 2 ]] ) )


## *********** fitted values *************
print( fitted( fitw2sls1e ) )
print( fitted( fitw2sls1e$eq[[ 1 ]] ) )

print( fitted( fitw2sls2 ) )
print( fitted( fitw2sls2$eq[[ 2 ]] ) )

print( fitted( fitw2sls3 ) )
print( fitted( fitw2sls3$eq[[ 1 ]] ) )

print( fitted( fitw2sls4e ) )
print( fitted( fitw2sls4e$eq[[ 2 ]] ) )

print( fitted( fitw2sls5 ) )
print( fitted( fitw2sls5$eq[[ 1 ]] ) )

print( fitted( fitw2slsd1 ) )
print( fitted( fitw2slsd1$eq[[ 2 ]] ) )

print( fitted( fitw2slsd2e ) )
print( fitted( fitw2slsd2e$eq[[ 1 ]] ) )

print( fitted( fitw2slsd3e ) )
print( fitted( fitw2slsd3e$eq[[ 2 ]] ) )


## *********** predicted values *************
predictData <- Kmenta
predictData$price <- Kmenta$price * 0.9
predictData$income <- Kmenta$income * 1.1

print( predict( fitw2sls1e, se.fit = TRUE, interval = "prediction" ) )
print( predict( fitw2sls1e$eq[[ 1 ]] ) )

print( predict( fitw2sls2, se.pred = TRUE, interval = "confidence",
   level = 0.999, data = predictData ) )
print( predict( fitw2sls2$eq[[ 2 ]] ) )

print( predict( fitw2sls3, se.pred = TRUE, interval = "prediction",
   level = 0.975 ) )
print( predict( fitw2sls3$eq[[ 1 ]] ) )

print( predict( fitw2sls4e, se.fit = TRUE, interval = "confidence",
   level = 0.25 ) )
print( predict( fitw2sls4e$eq[[ 2 ]] ) )

print( predict( fitw2sls5, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, data = predictData ) )
print( predict( fitw2sls5$eq[[ 1 ]] ) )

print( predict( fitw2slsd1, se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.99 ) )
print( predict( fitw2slsd1$eq[[ 2 ]] ) )

print( predict( fitw2slsd2e, se.fit = TRUE, interval = "prediction",
   level = 0.9, data = predictData ) )
print( predict( fitw2slsd2e$eq[[ 1 ]] ) )

print( predict( fitw2slsd3e, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.01 ) )
print( predict( fitw2slsd3e$eq[[ 2 ]] ) )

# predict just one observation
smallData <- data.frame( price = 130, income = 150, farmPrice = 120,
   trend = 25, consump = 1 ) ### consump should be removed later!!!!

print( predict( fitw2sls1e, data = smallData ) )
print( predict( fitw2sls1e$eq[[ 1 ]], data = smallData ) )

print( predict( fitw2sls2, se.fit = TRUE, level = 0.9,
   data = smallData ) )
print( predict( fitw2sls2$eq[[ 1 ]], se.pred = TRUE, level = 0.99,
   data = smallData ) )

print( predict( fitw2sls3, interval = "prediction", level = 0.975,
   data = smallData ) )
print( predict( fitw2sls3$eq[[ 1 ]], interval = "confidence", level = 0.8,
   data = smallData ) )

print( predict( fitw2sls4e, se.fit = TRUE, interval = "confidence",
   level = 0.999, data = smallData ) )
print( predict( fitw2sls4e$eq[[ 2 ]], se.pred = TRUE, interval = "prediction",
   level = 0.75, data = smallData ) )

print( predict( fitw2sls5, se.fit = TRUE, interval = "prediction",
   data = smallData ) )
print( predict( fitw2sls5$eq[[ 1 ]], se.pred = TRUE, interval = "confidence",
   data = smallData ) )

print( predict( fitw2slsd2e, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, data = smallData ) )
print( predict( fitw2slsd2e$eq[[ 1 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.25, data = smallData ) )


## ************ correlation of predicted values ***************
print( correlation.systemfit( fitw2sls1e, 1, 2 ) )

print( correlation.systemfit( fitw2sls2, 2, 1 ) )

print( correlation.systemfit( fitw2sls3, 1, 2 ) )

print( correlation.systemfit( fitw2sls4e, 2, 1 ) )

print( correlation.systemfit( fitw2sls5, 1, 2 ) )

print( correlation.systemfit( fitw2slsd1, 2, 1 ) )

print( correlation.systemfit( fitw2slsd2e, 1, 2 ) )

print( correlation.systemfit( fitw2slsd3e, 2, 1 ) )


## ************** F tests ****************
# testing first restriction
print( ftest.systemfit( fitw2sls1, restrm ) )
print( ftest.systemfit( fitw2slsd1e, restrm ) )

# testing second restriction
restrOnly2m <- matrix(0,1,7)
restrOnly2q <- 0.5
restrOnly2m[1,2] <- -1
restrOnly2m[1,5] <-  1
# first restriction not imposed 
print( ftest.systemfit( fitw2sls1e, restrOnly2m, restrOnly2q ) )
print( ftest.systemfit( fitw2slsd1, restrOnly2m, restrOnly2q ) )
# first restriction imposed
print( ftest.systemfit( fitw2sls2, restrOnly2m, restrOnly2q ) )
print( ftest.systemfit( fitw2sls3, restrOnly2m, restrOnly2q ) )
print( ftest.systemfit( fitw2slsd2e, restrOnly2m, restrOnly2q ) )
print( ftest.systemfit( fitw2slsd3e, restrOnly2m, restrOnly2q ) )

# testing both of the restrictions
print( ftest.systemfit( fitw2sls1e, restr2m, restr2q ) )
print( ftest.systemfit( fitw2slsd1, restr2m, restr2q ) )


## ************** Wald tests ****************
# testing first restriction
print( waldtest.systemfit( fitw2sls1, restrm ) )
print( waldtest.systemfit( fitw2slsd1e, restrm ) )

# testing second restriction
# first restriction not imposed
print( waldtest.systemfit( fitw2sls1e, restrOnly2m, restrOnly2q ) )
print( waldtest.systemfit( fitw2slsd1, restrOnly2m, restrOnly2q ) )
# first restriction imposed
print( waldtest.systemfit( fitw2sls2, restrOnly2m, restrOnly2q ) )
print( waldtest.systemfit( fitw2sls3, restrOnly2m, restrOnly2q ) )
print( waldtest.systemfit( fitw2slsd2e, restrOnly2m, restrOnly2q ) )
print( waldtest.systemfit( fitw2slsd3e, restrOnly2m, restrOnly2q ) )

# testing both of the restrictions
print( waldtest.systemfit( fitw2sls1e, restr2m, restr2q ) )
print( waldtest.systemfit( fitw2slsd1, restr2m, restr2q ) )

