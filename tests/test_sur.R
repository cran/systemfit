
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

# the standard equations do not converge and lead to a singular weighting matrix
# both in R and in EViews, since both equations have the same endogenous variable
supply2 <- price ~ income + farmPrice + trend
system2 <- list( demand, supply2 )


## *************** SUR estimation ************************
fitsur1 <- systemfit( "SUR", system, labels, data = Kmenta )
print( summary( fitsur1 ) )
print( round( fitsur1$bcov, digits = 6 ) )

## ********************* SUR (EViews-like) *****************
fitsur1e <- systemfit( "SUR", system, labels, data = Kmenta, rcovformula = 0,
   probdfsys = TRUE )
print( summary( fitsur1e ) )
print( round( fitsur1e$bcov, digits = 6 ) )

## ********************* SUR (rcovformula=2) *****************
fitsur1c <- systemfit( "SUR", system, labels, data = Kmenta, rcovformula = 2 )
print( summary( fitsur1c ) )
print( round( fitsur1c$bcov, digits = 6 ) )

## *************** SUR (rcovformula=2, probdfsys = TRUE ) ***************
fitsur1cp <- systemfit( "SUR", system, labels, data = Kmenta, rcovformula = 2,
   probdfsys = TRUE )
print( summary( fitsur1cp ) )
print( round( fitsur1cp$bcov, digits = 6 ) )

## ********************* SUR (rcovformula=3) *****************
fitsur1c <- systemfit( "SUR", system, labels, data = Kmenta, rcovformula = 3 )
print( summary( fitsur1c ) )
print( round( fitsur1c$bcov, digits = 6 ) )

## *************** SUR with cross-equation restriction **************
fitsur2 <- systemfit( "SUR", system, labels, data = Kmenta, R.restr = restrm )
print( summary( fitsur2 ) )
print( round( fitsur2$bcov, digits = 6 ) )

## *************** SUR with cross-equation restriction (EViews-like) **
fitsur2e <- systemfit( "SUR", system, labels, data = Kmenta, R.restr = restrm,
   rcovformula = 0 )
print( summary( fitsur2e ) )
print( round( fitsur2e$bcov, digits = 6 ) )

## *************** SUR with restriction via TX *******************
fitsur3 <- systemfit( "SUR", system, labels, data = Kmenta, TX = tc )
print( summary( fitsur3 ) )
print( round( fitsur3$bcov, digits = 6 ) )

## *************** SUR with restriction via TX (EViews-like) **************
fitsur3e <- systemfit( "SUR", system, labels, data = Kmenta, TX = tc,
   rcovformula = 0 )
print( summary( fitsur3e ) )
print( round( fitsur3e$bcov, digits = 6 ) )

## *************** SUR with 2 restrictions ***************************
fitsur4 <- systemfit( "SUR", system, labels, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q )
print( summary( fitsur4 ) )
print( round( fitsur4$bcov, digits = 6 ) )

## *************** SUR with 2 restrictions (EViews-like) **************
fitsur4e <- systemfit( "SUR", system, labels, data = Kmenta, rcovformula = 0,
   R.restr = restr2m, q.restr = restr2q )
print( summary( fitsur4e ) )
print( round( fitsur4e$bcov, digits = 6 ) )

## *************** SUR with 2 restrictions (rcovformula = 2) **************
fitsur4e <- systemfit( "SUR", system, labels, data = Kmenta, rcovformula = 2,
   R.restr = restr2m, q.restr = restr2q )
print( summary( fitsur4e ) )
print( round( fitsur4e$bcov, digits = 6 ) )

## *************** SUR with 2 restrictions (rcovformula = 3) **************
fitsur4e <- systemfit( "SUR", system, labels, data = Kmenta, rcovformula = 3,
   R.restr = restr2m, q.restr = restr2q )
print( summary( fitsur4e ) )
print( round( fitsur4e$bcov, digits = 6 ) )

## *************** SUR with 2 restrictions via R and TX ****************
fitsur5 <- systemfit( "SUR", system, labels, data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc )
print( summary( fitsur5 ) )
print( round( fitsur5$bcov, digits = 6 ) )

## *************** SUR with 2 restrictions via R and TX (EViews-like) **************
fitsur5e <- systemfit( "SUR", system, labels, data = Kmenta, rcovformula = 0,
   R.restr = restr3m, q.restr = restr3q, TX = tc )
print( summary( fitsur5e ) )
print( round( fitsur5e$bcov, digits = 6 ) )

## ************** iterated SUR ****************************
fitsuri1 <- systemfit( "SUR", system2, labels, data = Kmenta, maxit = 100 )
print( summary( fitsuri1 ) )
print( round( fitsuri1$bcov, digits = 6 ) )

## ************** iterated SUR (EViews-like) *****************
fitsuri1e <- systemfit( "SUR", system2, labels, data = Kmenta, rcovformula = 0,
   probdfsys = TRUE, maxit = 100 )
print( summary( fitsuri1e ) )
print( round( fitsuri1e$bcov, digits = 6 ) )

## ************** iterated SUR (rcovformula = 2) ****************************
fitsuri1c <- systemfit( "SUR", system2, labels, data = Kmenta, maxit = 100,
   rcovformula = 2 )
print( summary( fitsuri1c ) )
print( round( fitsuri1c$bcov, digits = 6 ) )

## ************** iterated SUR (rcovformula=2, probdfsys=TRUE) *****************
fitsuri1cp <- systemfit( "SUR", system2, labels, data = Kmenta, rcovformula = 2,
   probdfsys = TRUE, maxit = 100 )
print( summary( fitsuri1cp ) )
print( round( fitsuri1cp$bcov, digits = 6 ) )

## ************** iterated SUR (rcovformula = 3) ****************************
fitsuri1c <- systemfit( "SUR", system2, labels, data = Kmenta, maxit = 100,
   rcovformula = 3 )
print( summary( fitsuri1c ) )
print( round( fitsuri1c$bcov, digits = 6 ) )

## *********** iterated SUR with restriction *******************
fitsuri2 <- systemfit( "SUR", system2, labels, data = Kmenta, R.restr = restrm,
   maxit = 100 )
print( summary( fitsuri2 ) )
print( round( fitsuri2$bcov, digits = 6 ) )

## *********** iterated SUR with restriction (EViews-like) ***************
fitsuri2e <- systemfit( "SUR", system2, labels, data = Kmenta, R.restr = restrm,
   rcovformula = 0, maxit = 100 )
print( summary( fitsuri2e ) )
print( round( fitsuri2e$bcov, digits = 6 ) )

## *********** iterated SUR with restriction via TX ********************
fitsuri3 <- systemfit( "SUR", system2, labels, data = Kmenta, TX = tc,
   maxit = 100 )
print( summary( fitsuri3 ) )
print( round( fitsuri3$bcov, digits = 6 ) )

## *********** iterated SUR with restriction via TX (EViews-like) ***************
fitsuri3e <- systemfit( "SUR", system2, labels, data = Kmenta, TX = tc,
   rcovformula = 0, maxit = 100 )
print( summary( fitsuri3e ) )
print( round( fitsuri3e$bcov, digits = 6 ) )

## *************** iterated SUR with 2 restrictions ***************************
fitsuri4 <- systemfit( "SUR", system, labels, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, maxit = 100 )
print( summary( fitsuri4 ) )
print( round( fitsuri4$bcov, digits = 6 ) )

## *************** iterated SUR with 2 restrictions (EViews-like) **************
fitsuri4e <- systemfit( "SUR", system, labels, data = Kmenta, rcovformula = 0,
   R.restr = restr2m, q.restr = restr2q, maxit = 100 )
print( summary( fitsuri4e ) )
print( round( fitsuri4e$bcov, digits = 6 ) )

## *************** iterated SUR with 2 restrictions via R and TX ****************
fitsuri5 <- systemfit( "SUR", system, labels, data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, maxit = 100 )
print( summary( fitsuri5 ) )
print( round( fitsuri5$bcov, digits = 6 ) )

## ********* iterated SUR with 2 restrictions via R and TX (EViews-like) **********
fitsuri5e <- systemfit( "SUR", system, labels, data = Kmenta, rcovformula = 0,
   R.restr = restr3m, q.restr = restr3q, TX = tc, maxit = 100 )
print( summary( fitsuri5e ) )
print( round( fitsuri5e$bcov, digits = 6 ) )

## ********* iterated SUR with 2 restrictions via R and TX (rcovformula=2) **********
fitsuri5e <- systemfit( "SUR", system, labels, data = Kmenta, rcovformula = 2,
   R.restr = restr3m, q.restr = restr3q, TX = tc, maxit = 100 )
print( summary( fitsuri5e ) )
print( round( fitsuri5e$bcov, digits = 6 ) )

## ********* iterated SUR with 2 restrictions via R and TX (rcovformula=3) **********
# fitsuri5e <- systemfit( "SUR", system, labels, data = Kmenta, rcovformula = 3,
#    R.restr = restr3m, q.restr = restr3q, TX = tc, maxit = 100 )
# print( summary( fitsuri5e ) )
# print( round( fitsuri5e$bcov, digits = 6 ) )
# disabled, because the estimation does not converge
