
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


## *************** 3SLS estimation ************************
fit3sls <- list()
formulas <- c( "GLS", "IV", "Schmidt", "GMM", "EViews" )
for( i in seq( along = formulas ) ) {
   fit3sls[[ i ]] <- list()

   print( "***************************************************" )
   print( paste( "3SLS formula:", formulas[ i ] ) )
   print( "************* 3SLS *********************************" )
   fit3sls[[ i ]]$e1 <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = inst, formula3sls = formulas[ i ] )
   print( fit3sls[[ i ]]$e1 )
   print( round( fit3sls[[ i ]]$e1$bcov, digits = 6 ) )

   print( "********************* 3SLS EViews-like *****************" )
   fit3sls[[ i ]]$e1e <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = inst, rcovformula = 0, probdfsys = TRUE, formula3sls = formulas[ i ] )
   print( fit3sls[[ i ]]$e1e )
   print( round( fit3sls[[ i ]]$e1e$bcov, digits = 6 ) )

   print( "********************* 3SLS with rcovformula = 2 *****************" )
   fit3sls[[ i ]]$e1c <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = inst, rcovformula = 2, probdfsys = TRUE, formula3sls = formulas[ i ] )
   print( fit3sls[[ i ]]$e1c )
   print( round( fit3sls[[ i ]]$e1c$bcov, digits = 6 ) )

   print( "*************** 3SLS with restriction *****************" )
   fit3sls[[ i ]]$e2 <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = inst, R.restr = restrm, formula3sls = formulas[ i ] )
   print( fit3sls[[ i ]]$e2 )
   print( round( fit3sls[[ i ]]$e2$bcov, digits = 6 ) )

   print( "************** 3SLS with restriction (EViews-like) *****************" )
   fit3sls[[ i ]]$e2e <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = inst, rcovformula = 0, probdfsys = TRUE, R.restr = restrm,
      formula3sls = formulas[ i ] )
   print( fit3sls[[ i ]]$e2e )
   print( round( fit3sls[[ i ]]$e2e$bcov, digits = 6 ) )

   print( "*************** 3SLS with restriction via TX ********************" )
   fit3sls[[ i ]]$e3 <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = inst, TX = tc, formula3sls = formulas[ i ] )
   print( fit3sls[[ i ]]$e3 )
   print( round( fit3sls[[ i ]]$e3$bcov, digits = 6 ) )

   print( "*************** 3SLS with restriction via TX (EViews-like) *******" )
   fit3sls[[ i ]]$e3e <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = inst, rcovformula = 0, probdfsys = TRUE, TX = tc,
      formula3sls = formulas[ i ] )
   print( fit3sls[[ i ]]$e3e )
   print( round( fit3sls[[ i ]]$e3e$bcov, digits = 6 ) )

   print( "*************** 3SLS with 2 restrictions **********************" )
   fit3sls[[ i ]]$e4 <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = inst, R.restr = restr2m, q.restr = restr2q,
      formula3sls = formulas[ i ] )
   print( fit3sls[[ i ]]$e4 )
   print( round( fit3sls[[ i ]]$e4$bcov, digits = 6 ) )

   print( "*************** 3SLS with 2 restrictions (EViews-like) ************" )
   fit3sls[[ i ]]$e4e <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = inst, rcovformula = 0, probdfsys = TRUE, R.restr = restr2m,
      q.restr = restr2q, formula3sls = formulas[ i ] )
   print( fit3sls[[ i ]]$e4e )
   print( round( fit3sls[[ i ]]$e4e$bcov, digits = 6 ) )

   print( "*************** 3SLS with 2 restrictions via R and TX **********" )
   fit3sls[[ i ]]$e5 <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = inst, TX = tc, R.restr = restr3m, q.restr = restr3q,
      formula3sls = formulas[ i ] )
   print( fit3sls[[ i ]]$e5 )
   print( round( fit3sls[[ i ]]$e5$bcov, digits = 6 ) )

   print( "******** 3SLS with 2 restrictions via R and TX (EViews-like)*****" )
   fit3sls[[ i ]]$e5e <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = inst, TX = tc, rcovformula = 0, probdfsys = TRUE,
      R.restr = restr3m, q.restr = restr3q, formula3sls = formulas[ i ] )
   print( fit3sls[[ i ]]$e5e )
   print( round( fit3sls[[ i ]]$e5e$bcov, digits = 6 ) )
}

## ******************** iterated 3SLS **********************
fit3slsi <- list()
formulas <- c( "GLS", "IV", "Schmidt", "GMM", "EViews" )
for( i in seq( along = formulas ) ) {
   fit3slsi[[ i ]] <- list()

   print( "***************************************************" )
   print( paste( "3SLS formula:", formulas[ i ] ) )
   print( "************* 3SLS *********************************" )
   fit3slsi[[ i ]]$e1 <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = inst, formula3sls = formulas[ i ], maxiter = 100 )
   print( fit3slsi[[ i ]]$e1 )
   print( round( fit3slsi[[ i ]]$e1$bcov, digits = 6 ) )

   print( "********************* iterated 3SLS EViews-like ****************" )
   fit3slsi[[ i ]]$e1e <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = inst, rcovformula = 0, probdfsys = TRUE, formula3sls = formulas[ i ],
      maxiter = 100  )
   print( fit3slsi[[ i ]]$e1e )
   print( round( fit3slsi[[ i ]]$e1e$bcov, digits = 6 ) )

   print( "************** iterated 3SLS with rcovformula = 2 **************" )
   fit3slsi[[ i ]]$e1c <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = inst, rcovformula = 2, probdfsys = TRUE, formula3sls = formulas[ i ],
      maxiter = 100  )
   print( fit3slsi[[ i ]]$e1c )
   print( round( fit3slsi[[ i ]]$e1c$bcov, digits = 6 ) )

   print( "******* iterated 3SLS with restriction *****************" )
   fit3slsi[[ i ]]$e2 <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = inst, R.restr = restrm, formula3sls = formulas[ i ], maxiter = 100 )
   print( fit3slsi[[ i ]]$e2 )
   print( round( fit3slsi[[ i ]]$e2$bcov, digits = 6 ) )

   print( "********* iterated 3SLS with restriction (EViews-like) *********" )
   fit3slsi[[ i ]]$e2e <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = inst, rcovformula = 0, probdfsys = TRUE, R.restr = restrm,
      formula3sls = formulas[ i ], maxiter = 100 )
   print( fit3slsi[[ i ]]$e2e )
   print( round( fit3slsi[[ i ]]$e2e$bcov, digits = 6 ) )

   print( "********* iterated 3SLS with restriction via TX *****************" )
   fit3slsi[[ i ]]$e3 <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = inst, TX = tc, formula3sls = formulas[ i ], maxiter = 100 )
   print( fit3slsi[[ i ]]$e3 )
   print( round( fit3slsi[[ i ]]$e3$bcov, digits = 6 ) )

   print( "********* iterated 3SLS with restriction via TX (EViews-like) ***" )
   fit3slsi[[ i ]]$e3e <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = inst, rcovformula = 0, probdfsys = TRUE, TX = tc,
      formula3sls = formulas[ i ], maxiter = 100 )
   print( fit3slsi[[ i ]]$e3e )
   print( round( fit3slsi[[ i ]]$e3e$bcov, digits = 6 ) )

   print( "******** iterated 3SLS with 2 restrictions *********************" )
   fit3slsi[[ i ]]$e4 <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = inst, R.restr = restr2m, q.restr = restr2q,
      formula3sls = formulas[ i ], maxiter = 100 )
   print( fit3slsi[[ i ]]$e4 )
   print( round( fit3slsi[[ i ]]$e4$bcov, digits = 6 ) )

   print( "********* iterated 3SLS with 2 restrictions (EViews-like) *******" )
   fit3slsi[[ i ]]$e4e <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = inst, rcovformula = 0, probdfsys = TRUE, R.restr = restr2m,
      q.restr = restr2q, formula3sls = formulas[ i ], maxiter = 100 )
   print( fit3slsi[[ i ]]$e4e )
   print( round( fit3slsi[[ i ]]$e4e$bcov, digits = 6 ) )

   print( "******** iterated 3SLS with 2 restrictions via R and TX *********" )
   fit3slsi[[ i ]]$e5 <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = inst, TX = tc, R.restr = restr3m, q.restr = restr3q,
      formula3sls = formulas[ i ], maxiter = 100 )
   print( fit3slsi[[ i ]]$e5 )
   print( round( fit3slsi[[ i ]]$e5$bcov, digits = 6 ) )

   print( "*** iterated 3SLS with 2 restrictions via R and TX (EViews-like)**" )
   fit3slsi[[ i ]]$e5e <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = inst, TX = tc, rcovformula = 0, probdfsys = TRUE,
      R.restr = restr3m, q.restr = restr3q, formula3sls = formulas[ i ],
      maxiter = 100  )
   print( fit3slsi[[ i ]]$e5e )
   print( round( fit3slsi[[ i ]]$e5e$bcov, digits = 6 ) )
}

## **************** 3SLS with different instruments *************
fit3slsd <- list()
formulas <- c( "GLS", "IV", "Schmidt", "GMM", "EViews" )
for( i in seq( along = formulas ) ) {
   fit3slsd[[ i ]] <- list()

   print( "***************************************************" )
   print( paste( "3SLS formula:", formulas[ i ] ) )
   print( "************* 3SLS with different instruments **************" )
   fit3slsd[[ i ]]$e1 <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = instlist, formula3sls = formulas[ i ] )
   print( fit3slsd[[ i ]]$e1 )
   print( round( fit3slsd[[ i ]]$e1$bcov, digits = 6 ) )

   print( "******* 3SLS with different instruments (EViews-like) **********" )
   fit3slsd[[ i ]]$e1e <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = instlist, rcovformula = 0, probdfsys = TRUE, formula3sls = formulas[ i ] )
   print( fit3slsd[[ i ]]$e1e )
   print( round( fit3slsd[[ i ]]$e1e$bcov, digits = 6 ) )

   print( "**** 3SLS with different instruments and rcovformula = 2 ***" )
   fit3slsd[[ i ]]$e1c <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = instlist, rcovformula = 2, probdfsys = TRUE, formula3sls = formulas[ i ] )
   print( fit3slsd[[ i ]]$e1c )
   print( round( fit3slsd[[ i ]]$e1c$bcov, digits = 6 ) )

   print( "******* 3SLS with different instruments and restriction ********" )
   fit3slsd[[ i ]]$e2 <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = instlist, R.restr = restrm, formula3sls = formulas[ i ] )
   print( fit3slsd[[ i ]]$e2 )
   print( round( fit3slsd[[ i ]]$e2$bcov, digits = 6 ) )

   print( "** 3SLS with different instruments and restriction (EViews-like) *" )
   fit3slsd[[ i ]]$e2e <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = instlist, rcovformula = 0, probdfsys = TRUE, R.restr = restrm,
      formula3sls = formulas[ i ] )
   print( fit3slsd[[ i ]]$e2e )
   print( round( fit3slsd[[ i ]]$e2e$bcov, digits = 6 ) )

   print( "** 3SLS with different instruments and restriction via TX *******" )
   fit3slsd[[ i ]]$e3 <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = instlist, TX = tc, formula3sls = formulas[ i ] )
   print( fit3slsd[[ i ]]$e3 )
   print( round( fit3slsd[[ i ]]$e3$bcov, digits = 6 ) )

   print( "3SLS with different instruments with restriction via TX (EViews-like)" )
   fit3slsd[[ i ]]$e3e <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = instlist, rcovformula = 0, probdfsys = TRUE, TX = tc,
      formula3sls = formulas[ i ] )
   print( fit3slsd[[ i ]]$e3e )
   print( round( fit3slsd[[ i ]]$e3e$bcov, digits = 6 ) )

   print( "****** 3SLS with different instruments and 2 restrictions *********" )
   fit3slsd[[ i ]]$e4 <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = instlist, R.restr = restr2m, q.restr = restr2q,
      formula3sls = formulas[ i ] )
   print( fit3slsd[[ i ]]$e4 )
   print( round( fit3slsd[[ i ]]$e4$bcov, digits = 6 ) )

   print( "** 3SLS with different instruments and 2 restrictions (EViews-like) *" )
   fit3slsd[[ i ]]$e4e <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = instlist, rcovformula = 0, probdfsys = TRUE, R.restr = restr2m,
      q.restr = restr2q, formula3sls = formulas[ i ] )
   print( fit3slsd[[ i ]]$e4e )
   print( round( fit3slsd[[ i ]]$e4$bcov, digits = 6 ) )

   print( " 3SLS with different instruments with 2 restrictions via R and TX" )
   fit3slsd[[ i ]]$e5 <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = instlist, TX = tc, R.restr = restr3m, q.restr = restr3q,
      formula3sls = formulas[ i ] )
   print( fit3slsd[[ i ]]$e5 )
   print( round( fit3slsd[[ i ]]$e5$bcov, digits = 6 ) )

   print( "3SLS with diff. instruments and 2 restr. via R and TX (EViews-like)" )
   fit3slsd[[ i ]]$e5e <- systemfit( "3SLS", system, labels, data = kmenta,
      inst = instlist, TX = tc, rcovformula = 0, probdfsys = TRUE,
      R.restr = restr3m, q.restr = restr3q, formula3sls = formulas[ i ] )
   print( fit3slsd[[ i ]]$e5e )
   print( round( fit3slsd[[ i ]]$e5e$bcov, digits = 6 ) )
}
