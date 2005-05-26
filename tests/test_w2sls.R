
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


## ********************* W2SLS *****************
fitw2sls1 <- systemfit( "W2SLS", system, labels, data = kmenta, inst = inst )
print( fitw2sls1 )
print( round( fitw2sls1$bcov, digits = 6 ) )

## ********************* W2SLS (EViews-like) *****************
fitw2sls1e <- systemfit( "W2SLS", system, labels, data = kmenta, inst = inst,
   rcovformula = 0, probdfsys = TRUE )
print( fitw2sls1e )
print( round( fitw2sls1e$bcov, digits = 6 ) )

## ********************* W2SLS with restriction *******************
fitw2sls2 <- systemfit( "W2SLS", system, labels, data = kmenta, R.restr = restrm,
   inst = inst )
print( fitw2sls2 )
print( round( fitw2sls2$bcov, digits = 6 ) )

## ********************* W2SLS with restriction (EViews-like) **************
fitw2sls2e <- systemfit( "W2SLS", system, labels, data = kmenta, R.restr = restrm,
   inst = inst, rcovformula = 0, probdfsys = TRUE )
print( fitw2sls2e )
print( round( fitw2sls2e$bcov, digits = 6 ) )

## ********************* W2SLS with restriction via TX *******************
fitw2sls3 <- systemfit( "W2SLS", system, labels, data = kmenta, TX = tc, inst = inst )
print( fitw2sls3 )
print( round( fitw2sls3$bcov, digits = 6 ) )

## ********************* W2SLS with restriction via TX (EViews-like) **************
fitw2sls3e <- systemfit( "W2SLS", system, labels, data = kmenta, TX = tc, inst = inst,
   rcovformula = 0, probdfsys = TRUE )
print( fitw2sls3e )
print( round( fitw2sls3e$bcov, digits = 6 ) )

## ***************** W2SLS with 2 restrictions ********************
fitw2sls4 <- systemfit( "W2SLS", system, labels, data = kmenta, R.restr = restr2m,
   q.restr = restr2q, inst = inst )
print( fitw2sls4 )
print( round( fitw2sls4$bcov, digits = 6 ) )

## ***************** W2SLS with 2 restrictions (EViews-like) **************
fitw2sls4e <- systemfit( "W2SLS", system, labels, data = kmenta, R.restr = restr2m,
   q.restr = restr2q, inst = inst, rcovformula = 0, probdfsys = TRUE )
print( fitw2sls4e )
print( round( fitw2sls4e$bcov, digits = 6 ) )

## ***************** W2SLS with 2 restrictions via R and TX ******************
fitw2sls5 <- systemfit( "W2SLS", system, labels, data = kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, inst = inst )
print( fitw2sls5 )
print( round( fitw2sls5$bcov, digits = 6 ) )

## ***************** W2SLS with 2 restrictions via R and TX (EViews-like) **************
fitw2sls5e <- systemfit( "W2SLS", system, labels, data = kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, inst = inst, rcovformula = 0, probdfsys = TRUE )
print( fitw2sls5e )
print( round( fitw2sls5e$bcov, digits = 6 ) )

## ****** 2SLS estimation with different instruments **********************
fitw2slsd1 <- systemfit( "W2SLS", system, labels, data = kmenta, inst = instlist )
print( fitw2slsd1 )
print( round( fitw2slsd1$bcov, digits = 6 ) )

## ****** 2SLS estimation with different instruments (EViews-like)******************
fitw2slsd1e <- systemfit( "W2SLS", system, labels, data = kmenta, inst = instlist,
   rcovformula = 0, probdfsys = TRUE )
print( fitw2slsd1e )
print( round( fitw2slsd1e$bcov, digits = 6 ) )

## **** W2SLS estimation with different instruments and restriction ********
fitw2slsd2 <- systemfit( "W2SLS", system, labels, data = kmenta, R.restr = restrm,
   inst = instlist )
print( fitw2slsd2 )
print( round( fitw2slsd2$bcov, digits = 6 ) )

## **** W2SLS estimation with different instruments and restriction (EViews-like)*
fitw2slsd2e <- systemfit( "W2SLS", system, labels, data = kmenta, R.restr = restrm,
   inst = instlist, rcovformula = 0, probdfsys = TRUE )
print( fitw2slsd2e )
print( round( fitw2slsd2e$bcov, digits = 6 ) )

## ** W2SLS estimation with different instruments and restriction via TX ****
fitw2slsd3 <- systemfit( "W2SLS", system, labels, data = kmenta, TX = tc,
   inst = instlist)
print( fitw2slsd3 )
print( round( fitw2slsd3$bcov, digits = 6 ) )

## W2SLS estimation with different instruments and restriction via TX (EViews-like)
fitw2slsd3e <- systemfit( "W2SLS", system, labels, data = kmenta, TX = tc,
   inst = instlist, rcovformula = 0, probdfsys = TRUE )
print( fitw2slsd3e )
print( round( fitw2slsd3e$bcov, digits = 6 ) )
