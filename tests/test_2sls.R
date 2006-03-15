
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
print( round( fit2sls1$bcov, digits = 6 ) )

## *************** 2SLS estimation (single.eq.sigma=F)*******************
fit2sls1s <- systemfit( "2SLS", system, labels, data = Kmenta, inst = inst,
   single.eq.sigma = FALSE )
print( summary( fit2sls1s ) )
print( round( fit2sls1s$bcov, digits = 6 ) )

## ********************* 2SLS (probdfsys = TRUE) *****************
fit2sls1p <- systemfit( "2SLS", system, labels, data = Kmenta, inst = inst,
   probdfsys = TRUE )
print( summary( fit2sls1p ) )
print( round( fit2sls1p$bcov, digits = 6 ) )

## ********************* 2SLS (rcovformula = 0) *****************
fit2sls1r <- systemfit( "2SLS", system, labels, data = Kmenta, inst = inst,
   rcovformula = 0 )
print( summary( fit2sls1r ) )
print( round( fit2sls1r$bcov, digits = 6 ) )

## *************** 2SLS (rcovformula=0, single.eq.sigma=F) *************
fit2sls1rs <- systemfit( "2SLS", system, labels, data = Kmenta, inst = inst,
   rcovformula = 0, single.eq.sigma = FALSE )
print( summary( fit2sls1rs ) )
print( round( fit2sls1rs$bcov, digits = 6 ) )

## ********************* 2SLS with restriction ********************
## **************** 2SLS with restriction (default)********************
fit2sls2 <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restrm,
   inst = inst )
print( summary( fit2sls2 ) )
print( round( fit2sls2$bcov, digits = 6 ) )

## ************* 2SLS with restriction (single.eq.sigma=T) *****************
fit2sls2s <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restrm,
   inst = inst, single.eq.sigma = TRUE )
print( summary( fit2sls2s ) )
print( round( fit2sls2s$bcov, digits = 6 ) )

## ********************* 2SLS with restriction (probdfsys=T) **************
fit2sls2p <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restrm,
   inst = inst, probdfsys = TRUE )
print( summary( fit2sls2p ) )
print( round( fit2sls2p$bcov, digits = 6 ) )

## ********************* 2SLS with restriction (rcovformula = 0) **************
fit2sls2r <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restrm,
   inst = inst, rcovformula = 0 )
print( summary( fit2sls2r ) )
print( round( fit2sls2r$bcov, digits = 6 ) )

## ******** 2SLS with restriction (rcovformula=0, single.eq.sigma=TRUE) *********
fit2sls2rs <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restrm,
   inst = inst, rcovformula = 0, single.eq.sigma = TRUE )
print( summary( fit2sls2rs ) )
print( round( fit2sls2rs$bcov, digits = 6 ) )

## ********************* 2SLS with restriction via TX ******************
## *************** 2SLS with restriction via TX (default )***************
fit2sls3 <- systemfit( "2SLS", system, labels, data = Kmenta, TX = tc,
   inst = inst, rcovformula = 0, probdfsys = TRUE )
print( summary( fit2sls3 ) )
print( round( fit2sls3$bcov, digits = 6 ) )

## ********************* 2SLS with restriction via TX (EViews-like) *******
fit2sls3e <- systemfit( "2SLS", system, labels, data = Kmenta, TX = tc,
   inst = inst, rcovformula = 0, probdfsys = TRUE )
print( summary( fit2sls3e ) )
print( round( fit2sls3e$bcov, digits = 6 ) )

## ***************** 2SLS with 2 restrictions *******************
## ************** 2SLS with 2 restrictions (default) **************
fit2sls4 <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, inst = inst )
print( summary( fit2sls4 ) )
print( round( fit2sls4$bcov, digits = 6 ) )

## ************ 2SLS with 2 restrictions (single.eq.sigma=T) **************
fit2sls4s <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, inst = inst, single.eq.sigma = TRUE )
print( summary( fit2sls4s ) )
print( round( fit2sls4s$bcov, digits = 6 ) )

## ***************** 2SLS with 2 restrictions (probdfsys=T) **************
fit2sls4p <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, inst = inst, probdfsys = TRUE )
print( summary( fit2sls4p ) )
print( round( fit2sls4p$bcov, digits = 6 ) )

## ***************** 2SLS with 2 restrictions (rcovformula=0) **************
fit2sls4r <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, inst = inst, rcovformula = 0 )
print( summary( fit2sls4r ) )
print( round( fit2sls4r$bcov, digits = 6 ) )

## ***** 2SLS with 2 restrictions (rcovformula=0, single.eq.sigma=T) *******
fit2sls4rs <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, inst = inst, rcovformula = 0, single.eq.sigma = TRUE )
print( summary( fit2sls4rs ) )
print( round( fit2sls4rs$bcov, digits = 6 ) )

## ************* 2SLS with 2 restrictions via R and TX ******************
## ******** 2SLS with 2 restrictions via R and TX (default) *************
fit2sls5 <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, inst = inst )
print( summary( fit2sls5 ) )
print( round( fit2sls5$bcov, digits = 6 ) )

## ******* 2SLS with 2 restrictions via R and TX (single.eq.sigma=T) ******
fit2sls5s <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, inst = inst, single.eq.sigma = TRUE )
print( summary( fit2sls5s ) )
print( round( fit2sls5s$bcov, digits = 6 ) )

## ********** 2SLS with 2 restrictions via R and TX (probdfsys=T) *******
fit2sls5p <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, inst = inst, probdfsys = TRUE )
print( summary( fit2sls5p ) )
print( round( fit2sls5p$bcov, digits = 6 ) )

## ************* 2SLS with 2 restrictions via R and TX (rcovformula=0) *********
fit2sls5r <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, inst = inst, rcovformula = 0 )
print( summary( fit2sls5r ) )
print( round( fit2sls5r$bcov, digits = 6 ) )

## ** 2SLS with 2 restrictions via R and TX (rcovformula=0, single.eq.sigma=T) **
fit2sls5rs <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, inst = inst, rcovformula = 0, single.eq.sigma = TRUE )
print( summary( fit2sls5rs ) )
print( round( fit2sls5rs$bcov, digits = 6 ) )

## *********** 2SLS estimation with different instruments **************
## ******* 2SLS estimation with different instruments (default) *********
fit2slsd1 <- systemfit( "2SLS", system, labels, data = Kmenta, inst = instlist )
print( summary( fit2slsd1 ) )
print( round( fit2slsd1$bcov, digits = 6 ) )

## *********** 2SLS estimation with different instruments (single.eq.sigma=F)*****
fit2slsd1s <- systemfit( "2SLS", system, labels, data = Kmenta, inst = instlist,
   single.eq.sigma = FALSE )
print( summary( fit2slsd1s ) )
print( round( fit2slsd1s$bcov, digits = 6 ) )

## ********* 2SLS estimation with different instruments (probdfsys=T) *******
fit2slsd1p <- systemfit( "2SLS", system, labels, data = Kmenta, inst = instlist,
   probdfsys = TRUE )
print( summary( fit2slsd1p ) )
print( round( fit2slsd1p$bcov, digits = 6 ) )

## ********* 2SLS estimation with different instruments (rcovformula=0) ******
fit2slsd1r <- systemfit( "2SLS", system, labels, data = Kmenta, inst = instlist,
   rcovformula = 0 )
print( summary( fit2slsd1r ) )
print( round( fit2slsd1r$bcov, digits = 6 ) )

## 2SLS estimation with different instruments (rcovformula=0,single.eq.sigma=F)
fit2slsd1r <- systemfit( "2SLS", system, labels, data = Kmenta, inst = instlist,
   rcovformula = 0, single.eq.sigma = FALSE )
print( summary( fit2slsd1r ) )
print( round( fit2slsd1r$bcov, digits = 6 ) )

## **** 2SLS estimation with different instruments and restriction *******
## ** 2SLS estimation with different instruments and restriction (default) ****
fit2slsd2 <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restrm,
   inst = instlist )
print( summary( fit2slsd2 ) )
print( round( fit2slsd2$bcov, digits = 6 ) )

## 2SLS estimation with different instruments and restriction (single.eq.sigma=T)
fit2slsd2s <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restrm,
   inst = instlist, single.eq.sigma = TRUE )
print( summary( fit2slsd2s ) )
print( round( fit2slsd2s$bcov, digits = 6 ) )

## **** 2SLS estimation with different instruments and restriction (probdfsys=F)
fit2slsd2p <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restrm,
   inst = instlist, probdfsys = FALSE )
print( summary( fit2slsd2p ) )
print( round( fit2slsd2p$bcov, digits = 6 ) )

## **** 2SLS estimation with different instruments and restriction (rcovformula=0)
fit2slsd2r <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restrm,
   inst = instlist, rcovformula = 0 )
print( summary( fit2slsd2r ) )
print( round( fit2slsd2r$bcov, digits = 6 ) )

## 2SLS estimation with different instr. and restr. (rcovformula=0, single.eq.sigma=T)
fit2slsd2rs <- systemfit( "2SLS", system, labels, data = Kmenta, R.restr = restrm,
   inst = instlist, rcovformula = 0, single.eq.sigma = TRUE )
print( summary( fit2slsd2rs ) )
print( round( fit2slsd2rs$bcov, digits = 6 ) )

## **** 2SLS estimation with different instruments and restriction via TX *
## 2SLS estimation with different instruments and restriction via TX (default)
fit2slsd3 <- systemfit( "2SLS", system, labels, data = Kmenta, TX = tc,
   inst = instlist )
print( summary( fit2slsd3 ) )
print( round( fit2slsd3$bcov, digits = 6 ) )

## **** 2SLS estimation with different instr. and restr. via TX (rcovformula=0)
fit2slsd3r <- systemfit( "2SLS", system, labels, data = Kmenta, TX = tc,
   inst = instlist, rcovformula = 0 )
print( summary( fit2slsd3r ) )
print( round( fit2slsd3r$bcov, digits = 6 ) )
