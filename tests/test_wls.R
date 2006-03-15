
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


## *************** WLS estimation ************************
fitwls1 <- systemfit( "WLS", system, labels, data = Kmenta )
print( summary( fitwls1 ) )
print( round( fitwls1$bcov, digits = 6 ) )

## *************** WLS estimation (EViews-like) ************************
fitwls1e <- systemfit( "WLS", system, labels, data = Kmenta, rcovformula = 0,
   probdfsys = TRUE )
print( summary( fitwls1e ) )
print( round( fitwls1e$bcov, digits = 6 ) )

## ************** WLS with cross-equation restriction ***************
fitwls2 <- systemfit( "WLS", system, labels, data = Kmenta, R.restr = restrm )
print( summary( fitwls2 ) )
print( round( fitwls2$bcov, digits = 6 ) )

## ************** WLS with cross-equation restriction (EViews-like) *******
fitwls2e <- systemfit( "WLS", system, labels, data = Kmenta, R.restr = restrm,
   rcovformula = 0 )
print( summary( fitwls2e ) )
print( round( fitwls2e$bcov, digits = 6 ) )

## ******* WLS with cross-equation restriction via TX **********
fitwls3 <- systemfit("WLS", system, labels, data = Kmenta, TX = tc,)
print( summary( fitwls3 ) )
print( round( fitwls3$bcov, digits = 6 ) )

## ******* WLS with cross-equation restriction via TX (EViews-like) *****
fitwls3e <- systemfit("WLS", system, labels, data = Kmenta, TX = tc,
   rcovformula = 0 )
print( summary( fitwls3e ) )
print( round( fitwls3e$bcov, digits = 6 ) )

## ***** WLS with 2 cross-equation restrictions ***************
fitwls4 <- systemfit("WLS", system, labels, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q )
print( summary( fitwls4 ) )
print( round( fitwls4$bcov, digits = 6 ) )

## ***** WLS with 2 cross-equation restrictions (EViews-like) **********
fitwls4e <- systemfit("WLS", system, labels, data = Kmenta, rcovformula = 0,
   R.restr = restr2m, q.restr = restr2q )
print( summary( fitwls4e ) )
print( round( fitwls4e$bcov, digits = 6 ) )

## *********** WLS with 2 cross-equation restrictions via R and TX ******
fitwls5 <- systemfit( "WLS", system, labels, data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc )
print( summary( fitwls5 ) )
print( round( fitwls5$bcov, digits = 6 ) )

## *********** WLS with 2 cross-equation restrictions via R and TX (EViews-like)
fitwls5e <- systemfit( "WLS", system, labels, data = Kmenta, rcovformula = 0,
   R.restr = restr3m, q.restr = restr3q, TX = tc )
print( summary( fitwls5e ) )
print( round( fitwls5e$bcov, digits = 6 ) )

## *************** iterated WLS estimation *********************
fitwlsi1 <- systemfit( "WLS", system, labels, data = Kmenta, probdfsys = TRUE,
   maxit = 100 )
print( summary( fitwlsi1 ) )
print( round( fitwlsi1$bcov, digits = 6 ) )

## *************** iterated WLS estimation (EViews-like) ************
fitwlsi1e <- systemfit( "WLS", system, labels, data = Kmenta, rcovformula = 0,
   probdfsys = TRUE, maxit = 100 )
print( summary( fitwlsi1e ) )
print( round( fitwlsi1e$bcov, digits = 6 ) )

## ****** iterated WLS with cross-equation restriction ***************
fitwlsi2 <- systemfit( "WLS", system, labels, data = Kmenta, R.restr = restrm,
   maxit = 100 )
print( summary( fitwlsi2 ) )
print( round( fitwlsi2$bcov, digits = 6 ) )

## ****** iterated WLS with cross-equation restriction (EViews-like) ********
fitwlsi2e <- systemfit( "WLS", system, labels, data = Kmenta, R.restr = restrm,
   rcovformula = 0, maxit = 100 )
print( summary( fitwlsi2e ) )
print( round( fitwlsi2e$bcov, digits = 6 ) )

## ******* iterated WLS with cross-equation restriction via TX **********
fitwlsi3 <- systemfit( "WLS", system, labels, data = Kmenta, TX = tc,
   maxit = 100 )
print( summary( fitwlsi3 ) )
print( round( fitwlsi3$bcov, digits = 6 ) )

## ******* iterated WLS with cross-equation restriction via TX (EViews-like) ***
fitwlsi3e <- systemfit( "WLS", system, labels, data = Kmenta, TX = tc,
   rcovformula = 0, maxit = 100 )
print( summary( fitwlsi3e ) )
print( round( fitwlsi3e$bcov, digits = 6 ) )

## ******* iterated WLS with 2 cross-equation restrictions ***********
fitwlsi4 <- systemfit( "WLS", system, labels, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, maxit = 100 )
print( summary( fitwlsi4 ) )
print( round( fitwlsi4$bcov, digits = 6 ) )

## ******* iterated WLS with 2 cross-equation restrictions (EViews-like) *****
fitwlsi4e <- systemfit( "WLS", system, labels, data = Kmenta, rcovformula = 0,
   R.restr = restr2m, q.restr = restr2q, maxit = 100 )
print( summary( fitwlsi4e ) )
print( round( fitwlsi4e$bcov, digits = 6 ) )

## ***** iterated WLS with 2 cross-equation restrictions via R and TX ******
fitwlsi5 <- systemfit( "WLS", system, labels, data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, maxit = 100 )
print( summary( fitwlsi5 ) )
print( round( fitwlsi5$bcov, digits = 6 ) )

## *** iterated WLS with 2 cross-equation restrictions via R and TX (EViews-like)
fitwlsi5e <- systemfit( "WLS", system, labels, data = Kmenta, rcovformula = 0,
   R.restr = restr3m, q.restr = restr3q, TX = tc, maxit = 100 )
print( summary( fitwlsi5e ) )
print( round( fitwlsi5e$bcov, digits = 6 ) )

