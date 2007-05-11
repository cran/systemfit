
library( systemfit )
data( "Kmenta" )

eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
inst <- ~ income + farmPrice + trend
eqSystem <- list( demand = eqDemand, supply = eqSupply )
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


## ******************* unrestricted estimation *****************
## ******************** default estimation *********************
fit2sls1 <- systemfit( "2SLS", eqSystem, data = Kmenta, inst = inst )
fit3sls1 <- systemfit( "3SLS", eqSystem, data = Kmenta, inst = inst )
print( hausman.systemfit( fit2sls1, fit3sls1 ) )

## ************** 2SLS estimation with single.eq.sigma = FALSE *****************
fit2sls1s <- systemfit( "2SLS", eqSystem, data = Kmenta, inst = inst,
   single.eq.sigma = FALSE )
print( hausman.systemfit( fit2sls1s, fit3sls1 ) )

## ******************* estimations with rcovformula = 0 *****************
fit2sls1r <- systemfit( "2SLS", eqSystem, data = Kmenta, inst = inst,
   rcovformula = 0 )
fit3sls1r <- systemfit( "3SLS", eqSystem, data = Kmenta, inst = inst,
   rcovformula = 0 )
print( hausman.systemfit( fit2sls1r, fit3sls1r ) )


## ********************* estimation with restriction ********************
## *********************** default estimation ***********************
fit2sls2 <- systemfit( "2SLS", eqSystem, data = Kmenta, R.restr = restrm,
   inst = inst )
fit3sls2 <- systemfit( "3SLS", eqSystem, data = Kmenta, R.restr = restrm,
   inst = inst )
# print( hausman.systemfit( fit2sls2, fit3sls2 ) )

## ************* 2SLS estimation with single.eq.sigma = TRUE *****************
fit2sls2s <- systemfit( "2SLS", eqSystem, data = Kmenta, R.restr = restrm,
   inst = inst, single.eq.sigma = TRUE )
# print( hausman.systemfit( fit2sls2s, fit3sls2 ) )

## ********************* estimations with rcovformula = 0 **************
fit2sls2r <- systemfit( "2SLS", eqSystem, data = Kmenta, R.restr = restrm,
   inst = inst, rcovformula = 0 )
fit3sls2r <- systemfit( "3SLS", eqSystem, data = Kmenta, R.restr = restrm,
   inst = inst, rcovformula = 0 )
# print( hausman.systemfit( fit2sls2r, fit3sls2r ) )


## ****************** estimation with restriction via TX ******************
## ********************** default estimation ********************
fit2sls3 <- systemfit( "2SLS", eqSystem, data = Kmenta, TX = tc,
   inst = inst )
fit3sls3 <- systemfit( "3SLS", eqSystem, data = Kmenta, TX = tc,
   inst = inst )
print( hausman.systemfit( fit2sls3, fit3sls3 ) )

## ******************* estimations with rcovformula = 0 *******
fit2sls3r <- systemfit( "2SLS", eqSystem, data = Kmenta, TX = tc,
   inst = inst, rcovformula = 0 )
fit3sls3r <- systemfit( "3SLS", eqSystem, data = Kmenta, TX = tc,
   inst = inst, rcovformula = 0 )
print( hausman.systemfit( fit2sls3r, fit3sls3r ) )


## ***************** estimations with 2 restrictions *******************
## *********************** default estimations **************
fit2sls4 <- systemfit( "2SLS", eqSystem, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, inst = inst )
fit3sls4 <- systemfit( "3SLS", eqSystem, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, inst = inst )
# print( hausman.systemfit( fit2sls4, fit3sls4 ) )

## ***************** estimations with rcovformula = 0 **************
fit2sls4r <- systemfit( "2SLS", eqSystem, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, inst = inst, rcovformula = 0 )
fit3sls4r <- systemfit( "3SLS", eqSystem, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, inst = inst, rcovformula = 0 )
# print( hausman.systemfit( fit2sls4r, fit3sls4r ) )


## *********** estimations with 2 restrictions via R and TX ***************
## ***************** default estimations *******************
fit2sls5 <- systemfit( "2SLS", eqSystem, data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, inst = inst )
fit3sls5 <- systemfit( "3SLS", eqSystem, data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, inst = inst )
# print( hausman.systemfit( fit2sls5, fit3sls5 ) )

## ************* estimations with rcovformula = 0 *********
fit2sls5r <- systemfit( "2SLS", eqSystem, data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, inst = inst, rcovformula = 0 )
fit3sls5r <- systemfit( "3SLS", eqSystem, data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, inst = inst, rcovformula = 0 )
# print( hausman.systemfit( fit2sls5r, fit3sls5r ) )
