
library( systemfit )
data( kmenta )

demand <- q ~ p + d
supply <- q ~ p + f + a
inst   <- ~ d + f + a
inst1  <- ~ d + f
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


## *************** WLS estimation ************************
fitwls1 <- systemfit( "WLS", system, labels, data = kmenta )
print( fitwls1 )
print( fitwls1$bcov )

## *************** WLS estimation (EViews-like) ************************
fitwls1e <- systemfit( "WLS", system, labels, data = kmenta, rcovformula = 0,
   probdfsys = TRUE )
print( fitwls1e )
print( fitwls1e$bcov )

## ************** WLS with cross-equation restriction ***************
fitwls2 <- systemfit( "WLS", system, labels, data = kmenta, R.restr = restrm )
print( fitwls2 )
print( fitwls2$bcov )

## ************** WLS with cross-equation restriction (EViews-like) *******
fitwls2e <- systemfit( "WLS", system, labels, data = kmenta, R.restr = restrm,
   rcovformula = 0 )
print( fitwls2e )
print( fitwls2e$bcov )

## ******* WLS with cross-equation restriction via TX **********
fitwls3 <- systemfit("WLS", system, labels, data = kmenta, TX = tc,)
print( fitwls3 )
print( fitwls3$bcov )

## ******* WLS with cross-equation restriction via TX (EViews-like) *****
fitwls3e <- systemfit("WLS", system, labels, data = kmenta, TX = tc,
   rcovformula = 0 )
print( fitwls3e )
print( fitwls3e$bcov )

## ***** WLS with 2 cross-equation restrictions ***************
fitwls4 <- systemfit("WLS", system, labels, data = kmenta, R.restr = restr2m,
   q.restr = restr2q )
print( fitwls4 )
print( fitwls4$bcov )

## ***** WLS with 2 cross-equation restrictions (EViews-like) **********
fitwls4e <- systemfit("WLS", system, labels, data = kmenta, rcovformula = 0,
   R.restr = restr2m, q.restr = restr2q )
print( fitwls4e )
print( fitwls4e$bcov )

## *********** WLS with 2 cross-equation restrictions via R and TX ******
fitwls5 <- systemfit( "WLS", system, labels, data = kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc )
print( fitwls5 )
print( fitwls5$bcov )

## *********** WLS with 2 cross-equation restrictions via R and TX (EViews-like)
fitwls5e <- systemfit( "WLS", system, labels, data = kmenta, rcovformula = 0,
   R.restr = restr3m, q.restr = restr3q, TX = tc )
print( fitwls5e )
print( fitwls5e$bcov )

## *************** iterated WLS estimation *********************
fitwlsi1 <- systemfit( "WLS", system, labels, data = kmenta, probdfsys = TRUE,
   maxit = 100 )
print( fitwlsi1 )
print( fitwlsi1$bcov )

## *************** iterated WLS estimation (EViews-like) ************
fitwlsi1e <- systemfit( "WLS", system, labels, data = kmenta, rcovformula = 0,
   probdfsys = TRUE, maxit = 100 )
print( fitwlsi1e )
print( fitwlsi1e$bcov )

## ****** iterated WLS with cross-equation restriction ***************
fitwlsi2 <- systemfit( "WLS", system, labels, data = kmenta, R.restr = restrm,
   maxit = 100 )
print( fitwlsi2 )
print( fitwlsi2$bcov )

## ****** iterated WLS with cross-equation restriction (EViews-like) ********
fitwlsi2e <- systemfit( "WLS", system, labels, data = kmenta, R.restr = restrm,
   rcovformula = 0, maxit = 100 )
print( fitwlsi2e )
print( fitwlsi2e$bcov )

## ******* iterated WLS with cross-equation restriction via TX **********
fitwlsi3 <- systemfit( "WLS", system, labels, data = kmenta, TX = tc,
   maxit = 100 )
print( fitwlsi3 )
print( fitwlsi3$bcov )

## ******* iterated WLS with cross-equation restriction via TX (EViews-like) ***
fitwlsi3e <- systemfit( "WLS", system, labels, data = kmenta, TX = tc,
   rcovformula = 0, maxit = 100 )
print( fitwlsi3e )
print( fitwlsi3e$bcov )

## ******* iterated WLS with 2 cross-equation restrictions ***********
fitwlsi4 <- systemfit( "WLS", system, labels, data = kmenta, R.restr = restr2m,
   q.restr = restr2q, maxit = 100 )
print( fitwlsi4 )
print( fitwlsi4$bcov )

## ******* iterated WLS with 2 cross-equation restrictions (EViews-like) *****
fitwlsi4e <- systemfit( "WLS", system, labels, data = kmenta, rcovformula = 0,
   R.restr = restr2m, q.restr = restr2q, maxit = 100 )
print( fitwlsi4e )
print( fitwlsi4e$bcov )

## ***** iterated WLS with 2 cross-equation restrictions via R and TX ******
fitwlsi5 <- systemfit( "WLS", system, labels, data = kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, maxit = 100 )
print( fitwlsi5 )
print( fitwlsi5$bcov )

## *** iterated WLS with 2 cross-equation restrictions via R and TX (EViews-like)
fitwlsi5e <- systemfit( "WLS", system, labels, data = kmenta, rcovformula = 0,
   R.restr = restr3m, q.restr = restr3q, TX = tc, maxit = 100 )
print( fitwlsi5e )
print( fitwlsi5e$bcov )

