
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

# It is not possible to estimate OLS with systemfit
# exactly as EViews does, because EViews uses
# rcovformula == 1 for the coefficient covariance matrix and
# rcovformula == 0 for the residual covariance matrix, while
# systemfit uses always the same formulas for both calculations.

## *************** OLS estimation ************************
## ********** OLS estimation (default) ********************
fitols1 <- systemfit( "OLS", system, labels, data = kmenta )
print( fitols1 )
print( fitols1$bcov )

## ********** OLS estimation (no single.eq sigma=F) ******************
fitols1s <- systemfit( "OLS", system, labels, data = kmenta,
   single.eq.sigma = FALSE )
print( fitols1s )
print( fitols1s$bcov )

## ****************  OLS (probdfsys=T) ***********************
fitols1p <- systemfit( "OLS", system, labels, data = kmenta,
   probdfsys = TRUE )
print( fitols1p )
print( fitols1p$bcov )

## ****************  OLS (rcovformula=0) ***********************
fitols1r <- systemfit( "OLS", system, labels, data = kmenta,
   rcovformula = 0 )
print( fitols1r )
print( fitols1r$bcov )

## ********  OLS (rcovformula=0, single.eq.sigma=F) ***********
fitols1rs <- systemfit( "OLS", system, labels, data = kmenta,
   rcovformula = 0, single.eq.sigma = FALSE )
print( fitols1rs )
print( fitols1rs$bcov )

## ****************  OLS (rcovformula=2 ) ***********************
fitols1r <- systemfit( "OLS", system, labels, data = kmenta,
   rcovformula = 2 )
print( fitols1r )
print( fitols1r$bcov )

## ****************  OLS (rcovformula=3) ***********************
fitols1r <- systemfit( "OLS", system, labels, data = kmenta,
   rcovformula = 3 )
print( fitols1r )
print( fitols1r$bcov )

## ********  OLS (rcovformula=3, single.eq.sigma=F) ***********
fitols1rs <- systemfit( "OLS", system, labels, data = kmenta,
   rcovformula = 3, single.eq.sigma = FALSE )
print( fitols1rs )
print( fitols1rs$bcov )


## ********* OLS with cross-equation restriction ************
## ****** OLS with cross-equation restriction (default) *********
fitols2 <- systemfit( "OLS", system, labels, data = kmenta,
   R.restr = restrm )
print( fitols2 )
print( fitols2$bcov )

## ****** OLS with cross-equation restriction (single.eq.sigma=T) *******
fitols2s <- systemfit( "OLS", system, labels, data = kmenta,
   R.restr = restrm, single.eq.sigma = TRUE )
print( fitols2s )
print( fitols2s$bcov )

## ****** OLS with cross-equation restriction (probdfsys=F) *******
fitols2p <- systemfit( "OLS", system, labels, data = kmenta,
   R.restr = restrm, probdfsys = FALSE )
print( fitols2p )
print( fitols2p$bcov )

## ****** OLS with cross-equation restriction (rcovformula=0) *******
fitols2r <- systemfit( "OLS", system, labels, data = kmenta,
   R.restr = restrm, rcovformula = 0 )
print( fitols2r )
print( fitols2r$bcov )

## ** OLS with cross-equation restriction (rcovformula=0,single.eq.sigma=T) ***
fitols2rs <- systemfit( "OLS", system, labels, data = kmenta,
   R.restr = restrm, rcovformula = 0 )
print( fitols2rs )
print( fitols2rs$bcov )

## *** OLS with cross-equation restriction via TX ***
## *** OLS with cross-equation restriction via TX (default) ***
fitols3 <- systemfit( "OLS", system, labels, data = kmenta, TX = tc )
print( fitols3 )
print( fitols3$bcov )

## *** OLS with cross-equation restriction via TX (single.eq.sigma=T) ***
fitols3s <- systemfit( "OLS", system, labels, data = kmenta,
   TX = tc, single.eq.sigma = TRUE )
print( fitols3s )
print( fitols3s$bcov )

## *** OLS with cross-equation restriction via TX (probdfsys=F) ***
fitols3p <- systemfit( "OLS", system, labels, data = kmenta,
   TX = tc, probdfsys = FALSE )
print( fitols3p )
print( fitols3p$bcov )

## *** OLS with cross-equation restriction via TX (rcovformula=0) ***
fitols3r <- systemfit( "OLS", system, labels, data = kmenta,
   TX = tc, rcovformula = 0 )
print( fitols3r )
print( fitols3r$bcov )

## OLS with cross-equation restriction via TX (rcovformula=0,single.eq.sigma=T)
fitols3rs <- systemfit( "OLS", system, labels, data = kmenta,
   TX = tc, rcovformula = 0, single.eq.sigma = TRUE )
print( fitols3rs )
print( fitols3rs$bcov )

## ********* OLS with 2 cross-equation restrictions ***********
## ********* OLS with 2 cross-equation restrictions (default) ***********
fitols4 <- systemfit( "OLS", system, labels, data = kmenta, R.restr = restr2m,
   q.restr = restr2q )
print( fitols4 )
print( fitols4$bcov )

## ****** OLS with 2 cross-equation restrictions (single.eq.sigma=T) *******
fitols4s <- systemfit( "OLS", system, labels, data = kmenta, R.restr = restr2m,
   q.restr = restr2q, single.eq.sigma = T )
print( fitols4s )
print( fitols4s$bcov )

## ****** OLS with 2 cross-equation restrictions (probdfsys=F) *******
fitols4p <- systemfit( "OLS", system, labels, data = kmenta, R.restr = restr2m,
   q.restr = restr2q, probdfsys = F )
print( fitols4p )
print( fitols4p$bcov )

## ****** OLS with 2 cross-equation restrictions (rcovformula=0) *******
fitols4r <- systemfit( "OLS", system, labels, data = kmenta, R.restr = restr2m,
   q.restr = restr2q, rcovformula = 0 )
print( fitols4r )
print( fitols4r$bcov )

## OLS with 2 cross-equation restrictions (rcovformula=0, single.eq.sigma=T) *
fitols4rs <- systemfit( "OLS", system, labels, data = kmenta, R.restr = restr2m,
   q.restr = restr2q, rcovformula = 0, single.eq.sigma = T )
print( fitols4rs )
print( fitols4rs$bcov )

## ***** OLS with 2 cross-equation restrictions via R and TX ****
## ***** OLS with 2 cross-equation restrictions via R and TX (default) ****
fitols5 <- systemfit( "OLS", system, labels, data = kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, rcovformula = 0)
print( fitols5 )
print( fitols5$bcov )

## ***** OLS with 2 cross-equation restrictions via R and TX (single.eq.sigma=T) ****
fitols5s <- systemfit( "OLS", system, labels, data = kmenta,R.restr = restr3m,
   q.restr = restr3q, TX = tc, single.eq.sigma = T )
print( fitols5s )
print( fitols5s$bcov )

## ***** OLS with 2 cross-equation restrictions via R and TX (probdfsys=F) ****
fitols5p <- systemfit( "OLS", system, labels, data = kmenta,R.restr = restr3m,
   q.restr = restr3q, TX = tc, probdfsys = F )
print( fitols5p )
print( fitols5p$bcov )

## ***** OLS with 2 cross-equation restrictions via R and TX (rcovformula=0) ****
fitols5r <- systemfit( "OLS", system, labels, data = kmenta,R.restr = restr3m,
   q.restr = restr3q, TX = tc, rcovformula = 0 )
print( fitols5r )
print( fitols5r$bcov )

## OLS with 2 cross-equation restr. via R and TX (rcovformula=0,single.eq.sigma=T)
fitols5rs <- systemfit( "OLS", system, labels, data = kmenta,R.restr = restr3m,
   q.restr = restr3q, TX = tc, rcovformula = 0, single.eq.sigma = T )
print( fitols5rs )
print( fitols5rs$bcov )
