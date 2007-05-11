
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


## *************** 3SLS estimation ************************
fit3sls <- list()
formulas <- c( "GLS", "IV", "Schmidt", "GMM", "EViews" )
for( i in seq( along = formulas ) ) {
   fit3sls[[ i ]] <- list()

   print( "***************************************************" )
   print( paste( "3SLS formula:", formulas[ i ] ) )
   print( "************* 3SLS *********************************" )
   fit3sls[[ i ]]$e1 <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = inst, formula3sls = formulas[ i ] )
   print( summary( fit3sls[[ i ]]$e1 ) )

   print( "********************* 3SLS EViews-like *****************" )
   fit3sls[[ i ]]$e1e <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = inst, rcovformula = 0, probdfsys = TRUE, formula3sls = formulas[ i ] )
   print( summary( fit3sls[[ i ]]$e1e ) )

   print( "********************* 3SLS with rcovformula = 2 *****************" )
   fit3sls[[ i ]]$e1c <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = inst, rcovformula = 2, probdfsys = TRUE, formula3sls = formulas[ i ] )
   print( summary( fit3sls[[ i ]]$e1c ) )

   print( "*************** 3SLS with restriction *****************" )
   fit3sls[[ i ]]$e2 <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = inst, R.restr = restrm, formula3sls = formulas[ i ] )
   print( summary( fit3sls[[ i ]]$e2 ) )

   print( "************** 3SLS with restriction (EViews-like) *****************" )
   fit3sls[[ i ]]$e2e <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = inst, rcovformula = 0, probdfsys = TRUE, R.restr = restrm,
      formula3sls = formulas[ i ] )
   print( summary( fit3sls[[ i ]]$e2e ) )

   print( "*************** 3SLS with restriction via TX ********************" )
   fit3sls[[ i ]]$e3 <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = inst, TX = tc, formula3sls = formulas[ i ] )
   print( summary( fit3sls[[ i ]]$e3 ) )

   print( "*************** 3SLS with restriction via TX (EViews-like) *******" )
   fit3sls[[ i ]]$e3e <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = inst, rcovformula = 0, probdfsys = TRUE, TX = tc,
      formula3sls = formulas[ i ] )
   print( summary( fit3sls[[ i ]]$e3e ) )

   print( "*************** 3SLS with 2 restrictions **********************" )
   fit3sls[[ i ]]$e4 <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = inst, R.restr = restr2m, q.restr = restr2q,
      formula3sls = formulas[ i ] )
   print( summary( fit3sls[[ i ]]$e4 ) )

   print( "*************** 3SLS with 2 restrictions (EViews-like) ************" )
   fit3sls[[ i ]]$e4e <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = inst, rcovformula = 0, probdfsys = TRUE, R.restr = restr2m,
      q.restr = restr2q, formula3sls = formulas[ i ] )
   print( summary( fit3sls[[ i ]]$e4e ) )

   print( "*************** 3SLS with 2 restrictions via R and TX **********" )
   fit3sls[[ i ]]$e5 <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = inst, TX = tc, R.restr = restr3m, q.restr = restr3q,
      formula3sls = formulas[ i ] )
   print( summary( fit3sls[[ i ]]$e5 ) )

   print( "******** 3SLS with 2 restrictions via R and TX (EViews-like)*****" )
   fit3sls[[ i ]]$e5e <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = inst, TX = tc, rcovformula = 0, probdfsys = TRUE,
      R.restr = restr3m, q.restr = restr3q, formula3sls = formulas[ i ] )
   print( summary( fit3sls[[ i ]]$e5e ) )
}

## ******************** iterated 3SLS **********************
fit3slsi <- list()
formulas <- c( "GLS", "IV", "Schmidt", "GMM", "EViews" )
for( i in seq( along = formulas ) ) {
   fit3slsi[[ i ]] <- list()

   print( "***************************************************" )
   print( paste( "3SLS formula:", formulas[ i ] ) )
   print( "************* 3SLS *********************************" )
   fit3slsi[[ i ]]$e1 <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = inst, formula3sls = formulas[ i ], maxiter = 100 )
   print( summary( fit3slsi[[ i ]]$e1 ) )

   print( "********************* iterated 3SLS EViews-like ****************" )
   fit3slsi[[ i ]]$e1e <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = inst, rcovformula = 0, probdfsys = TRUE, formula3sls = formulas[ i ],
      maxiter = 100  )
   print( summary( fit3slsi[[ i ]]$e1e ) )

   print( "************** iterated 3SLS with rcovformula = 2 **************" )
   fit3slsi[[ i ]]$e1c <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = inst, rcovformula = 2, probdfsys = TRUE, formula3sls = formulas[ i ],
      maxiter = 100  )
   print( summary( fit3slsi[[ i ]]$e1c ) )

   print( "******* iterated 3SLS with restriction *****************" )
   fit3slsi[[ i ]]$e2 <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = inst, R.restr = restrm, formula3sls = formulas[ i ], maxiter = 100 )
   print( summary( fit3slsi[[ i ]]$e2 ) )

   print( "********* iterated 3SLS with restriction (EViews-like) *********" )
   fit3slsi[[ i ]]$e2e <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = inst, rcovformula = 0, probdfsys = TRUE, R.restr = restrm,
      formula3sls = formulas[ i ], maxiter = 100 )
   print( summary( fit3slsi[[ i ]]$e2e ) )

   print( "********* iterated 3SLS with restriction via TX *****************" )
   fit3slsi[[ i ]]$e3 <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = inst, TX = tc, formula3sls = formulas[ i ], maxiter = 100 )
   print( summary( fit3slsi[[ i ]]$e3 ) )

   print( "********* iterated 3SLS with restriction via TX (EViews-like) ***" )
   fit3slsi[[ i ]]$e3e <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = inst, rcovformula = 0, probdfsys = TRUE, TX = tc,
      formula3sls = formulas[ i ], maxiter = 100 )
   print( summary( fit3slsi[[ i ]]$e3e ) )

   print( "******** iterated 3SLS with 2 restrictions *********************" )
   fit3slsi[[ i ]]$e4 <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = inst, R.restr = restr2m, q.restr = restr2q,
      formula3sls = formulas[ i ], maxiter = 100 )
   print( summary( fit3slsi[[ i ]]$e4 ) )

   print( "********* iterated 3SLS with 2 restrictions (EViews-like) *******" )
   fit3slsi[[ i ]]$e4e <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = inst, rcovformula = 0, probdfsys = TRUE, R.restr = restr2m,
      q.restr = restr2q, formula3sls = formulas[ i ], maxiter = 100 )
   print( summary( fit3slsi[[ i ]]$e4e ) )

   print( "******** iterated 3SLS with 2 restrictions via R and TX *********" )
   fit3slsi[[ i ]]$e5 <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = inst, TX = tc, R.restr = restr3m, q.restr = restr3q,
      formula3sls = formulas[ i ], maxiter = 100 )
   print( summary( fit3slsi[[ i ]]$e5 ) )

   print( "*** iterated 3SLS with 2 restrictions via R and TX (EViews-like)**" )
   fit3slsi[[ i ]]$e5e <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = inst, TX = tc, rcovformula = 0, probdfsys = TRUE,
      R.restr = restr3m, q.restr = restr3q, formula3sls = formulas[ i ],
      maxiter = 100  )
   print( summary( fit3slsi[[ i ]]$e5e ) )
}

## **************** 3SLS with different instruments *************
fit3slsd <- list()
formulas <- c( "GLS", "IV", "Schmidt", "GMM", "EViews" )
for( i in seq( along = formulas ) ) {
   fit3slsd[[ i ]] <- list()

   print( "***************************************************" )
   print( paste( "3SLS formula:", formulas[ i ] ) )
   print( "************* 3SLS with different instruments **************" )
   fit3slsd[[ i ]]$e1 <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = instlist, formula3sls = formulas[ i ] )
   print( summary( fit3slsd[[ i ]]$e1 ) )

   print( "******* 3SLS with different instruments (EViews-like) **********" )
   fit3slsd[[ i ]]$e1e <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = instlist, rcovformula = 0, probdfsys = TRUE, formula3sls = formulas[ i ] )
   print( summary( fit3slsd[[ i ]]$e1e ) )

   print( "**** 3SLS with different instruments and rcovformula = 2 ***" )
   fit3slsd[[ i ]]$e1c <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = instlist, rcovformula = 2, probdfsys = TRUE, formula3sls = formulas[ i ] )
   print( summary( fit3slsd[[ i ]]$e1c ) )

   print( "******* 3SLS with different instruments and restriction ********" )
   fit3slsd[[ i ]]$e2 <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = instlist, R.restr = restrm, formula3sls = formulas[ i ] )
   print( summary( fit3slsd[[ i ]]$e2 ) )

   print( "** 3SLS with different instruments and restriction (EViews-like) *" )
   fit3slsd[[ i ]]$e2e <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = instlist, rcovformula = 0, probdfsys = TRUE, R.restr = restrm,
      formula3sls = formulas[ i ] )
   print( summary( fit3slsd[[ i ]]$e2e ) )

   print( "** 3SLS with different instruments and restriction via TX *******" )
   fit3slsd[[ i ]]$e3 <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = instlist, TX = tc, formula3sls = formulas[ i ] )
   print( summary( fit3slsd[[ i ]]$e3 ) )

   print( "3SLS with different instruments with restriction via TX (EViews-like)" )
   fit3slsd[[ i ]]$e3e <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = instlist, rcovformula = 0, probdfsys = TRUE, TX = tc,
      formula3sls = formulas[ i ] )
   print( summary( fit3slsd[[ i ]]$e3e ) )

   print( "****** 3SLS with different instruments and 2 restrictions *********" )
   fit3slsd[[ i ]]$e4 <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = instlist, R.restr = restr2m, q.restr = restr2q,
      formula3sls = formulas[ i ] )
   print( summary( fit3slsd[[ i ]]$e4 ) )

   print( "** 3SLS with different instruments and 2 restrictions (EViews-like) *" )
   fit3slsd[[ i ]]$e4e <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = instlist, rcovformula = 0, probdfsys = TRUE, R.restr = restr2m,
      q.restr = restr2q, formula3sls = formulas[ i ] )
   print( summary( fit3slsd[[ i ]]$e4e ) )

   print( " 3SLS with different instruments with 2 restrictions via R and TX" )
   fit3slsd[[ i ]]$e5 <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = instlist, TX = tc, R.restr = restr3m, q.restr = restr3q,
      formula3sls = formulas[ i ] )
   print( summary( fit3slsd[[ i ]]$e5 ) )

   print( "3SLS with diff. instruments and 2 restr. via R and TX (EViews-like)" )
   fit3slsd[[ i ]]$e5e <- systemfit( "3SLS", system, labels, data = Kmenta,
      inst = instlist, TX = tc, rcovformula = 0, probdfsys = TRUE,
      R.restr = restr3m, q.restr = restr3q, formula3sls = formulas[ i ] )
   print( summary( fit3slsd[[ i ]]$e5e ) )
}


## ****************** residuals **************************
print( residuals( fit3sls[[ 1 ]]$e1c ) )
print( residuals( fit3sls[[ 1 ]]$e1c$eq[[ 1 ]] ) )

print( residuals( fit3sls[[ 2 ]]$e2e ) )
print( residuals( fit3sls[[ 2 ]]$e2e$eq[[ 2 ]] ) )

print( residuals( fit3sls[[ 3 ]]$e3 ) )
print( residuals( fit3sls[[ 3 ]]$e3$eq[[ 1 ]] ) )

print( residuals( fit3sls[[ 4 ]]$e4e ) )
print( residuals( fit3sls[[ 4 ]]$e4e$eq[[ 2 ]] ) )

print( residuals( fit3sls[[ 5 ]]$e5 ) )
print( residuals( fit3sls[[ 5 ]]$e5$eq[[ 1 ]] ) )

print( residuals( fit3slsi[[ 2 ]]$e3e ) )
print( residuals( fit3slsi[[ 2 ]]$e3e$eq[[ 1 ]] ) )

print( residuals( fit3slsd[[ 3 ]]$e4 ) )
print( residuals( fit3slsd[[ 3 ]]$e4$eq[[ 2 ]] ) )


## *********** variance covariance matrix of the coefficients *******
print( round( vcov( fit3sls[[ 3 ]]$e1c ), digits = 6 ) )
print( round( vcov( fit3sls[[ 4 ]]$e1c$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fit3sls[[ 4 ]]$e2 ), digits = 6 ) )
print( round( vcov( fit3sls[[ 5 ]]$e2$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fit3sls[[ 5 ]]$e3e ), digits = 6 ) )
print( round( vcov( fit3sls[[ 1 ]]$e3e$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fit3sls[[ 1 ]]$e4 ), digits = 6 ) )
print( round( vcov( fit3sls[[ 2 ]]$e4$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fit3sls[[ 2 ]]$e5e ), digits = 6 ) )
print( round( vcov( fit3sls[[ 3 ]]$e5e$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fit3slsi[[ 4 ]]$e1e ), digits = 6 ) )
print( round( vcov( fit3slsi[[ 3 ]]$e1e$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fit3slsi[[ 5 ]]$e2e ), digits = 6 ) )
print( round( vcov( fit3slsi[[ 4 ]]$e2e$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fit3slsi[[ 1 ]]$e3 ), digits = 6 ) )
print( round( vcov( fit3slsi[[ 5 ]]$e3$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fit3slsi[[ 2 ]]$e4e ), digits = 6 ) )
print( round( vcov( fit3slsi[[ 1 ]]$e4e$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fit3slsi[[ 3 ]]$e5 ), digits = 6 ) )
print( round( vcov( fit3slsi[[ 2 ]]$e5$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fit3slsd[[ 5 ]]$e1c ), digits = 6 ) )
print( round( vcov( fit3slsd[[ 2 ]]$e1c$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fit3slsd[[ 1 ]]$e2 ), digits = 6 ) )
print( round( vcov( fit3slsd[[ 3 ]]$e2$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fit3slsd[[ 2 ]]$e3 ), digits = 6 ) )
print( round( vcov( fit3slsd[[ 4 ]]$e3$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fit3slsd[[ 3 ]]$e4 ), digits = 6 ) )
print( round( vcov( fit3slsd[[ 5 ]]$e4$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fit3slsd[[ 4 ]]$e5e ), digits = 6 ) )
print( round( vcov( fit3slsd[[ 1 ]]$e5e$eq[[ 2 ]] ), digits = 6 ) )


## *********** confidence intervals of coefficients *************
print( confint( fit3sls[[ 1 ]]$e1c ) )
print( confint( fit3sls[[ 1 ]]$e1c$eq[[ 1 ]], level = 0.9 ) )

print( confint( fit3sls[[ 2 ]]$e2e, level = 0.9 ) )
print( confint( fit3sls[[ 2 ]]$e2e$eq[[ 2 ]], level = 0.99 ) )

print( confint( fit3sls[[ 3 ]]$e3, level = 0.99 ) )
print( confint( fit3sls[[ 3 ]]$e3$eq[[ 1 ]], level = 0.5 ) )

print( confint( fit3sls[[ 4 ]]$e4e, level = 0.5 ) )
print( confint( fit3sls[[ 4 ]]$e4e$eq[[ 2 ]], level = 0.25 ) )

print( confint( fit3sls[[ 5 ]]$e5, level = 0.25 ) )
print( confint( fit3sls[[ 5 ]]$e5$eq[[ 1 ]], level = 0.975 ) )

print( confint( fit3slsi[[ 2 ]]$e3e, level = 0.975 ) )
print( confint( fit3slsi[[ 2 ]]$e3e$eq[[ 1 ]], level = 0.999 ) )

print( confint( fit3slsd[[ 3 ]]$e4, level = 0.999 ) )
print( confint( fit3slsd[[ 3 ]]$e4$eq[[ 2 ]] ) )


## *********** fitted values *************
print( fitted( fit3sls[[ 2 ]]$e1c ) )
print( fitted( fit3sls[[ 2 ]]$e1c$eq[[ 1 ]] ) )

print( fitted( fit3sls[[ 3 ]]$e2e ) )
print( fitted( fit3sls[[ 3 ]]$e2e$eq[[ 2 ]] ) )

print( fitted( fit3sls[[ 4 ]]$e3 ) )
print( fitted( fit3sls[[ 4 ]]$e3$eq[[ 1 ]] ) )

print( fitted( fit3sls[[ 5 ]]$e4e ) )
print( fitted( fit3sls[[ 5 ]]$e4e$eq[[ 2 ]] ) )

print( fitted( fit3sls[[ 1 ]]$e5 ) )
print( fitted( fit3sls[[ 1 ]]$e5$eq[[ 1 ]] ) )

print( fitted( fit3slsi[[ 3 ]]$e3e ) )
print( fitted( fit3slsi[[ 3 ]]$e3e$eq[[ 1 ]] ) )

print( fitted( fit3slsd[[ 4 ]]$e4 ) )
print( fitted( fit3slsd[[ 4 ]]$e4$eq[[ 2 ]] ) )


## *********** predicted values *************
predictData <- Kmenta
predictData$price <- Kmenta$price * 0.9
predictData$income <- Kmenta$income * 1.1

print( predict( fit3sls[[ 2 ]]$e1c, se.fit = TRUE, interval = "prediction" ) )
print( predict( fit3sls[[ 2 ]]$e1c$eq[[ 1 ]] ) )

print( predict( fit3sls[[ 3 ]]$e2e, se.pred = TRUE, interval = "confidence",
   level = 0.999, data = predictData ) )
print( predict( fit3sls[[ 3 ]]$e2e$eq[[ 2 ]] ) )

print( predict( fit3sls[[ 4 ]]$e3, se.pred = TRUE, interval = "prediction",
   level = 0.975 ) )
print( predict( fit3sls[[ 4 ]]$e3$eq[[ 1 ]] ) )

print( predict( fit3sls[[ 5 ]]$e4e, se.fit = TRUE, interval = "confidence",
   level = 0.25 ) )
print( predict( fit3sls[[ 5 ]]$e4e$eq[[ 2 ]] ) )

print( predict( fit3sls[[ 1 ]]$e5, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, data = predictData ) )
print( predict( fit3sls[[ 1 ]]$e5$eq[[ 1 ]] ) )

print( predict( fit3slsi[[ 3 ]]$e3e, se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.99 ) )
print( predict( fit3slsi[[ 3 ]]$e3e$eq[[ 1 ]] ) )

print( predict( fit3slsd[[ 4 ]]$e4, se.fit = TRUE, interval = "prediction",
   level = 0.9, data = predictData ) )
print( predict( fit3slsd[[ 4 ]]$e4$eq[[ 2 ]] ) )

# predict just one observation
smallData <- data.frame( price = 130, income = 150, farmPrice = 120,
   trend = 25, consump = 1 ) ### consump should be removed later!!!!

print( predict( fit3sls[[ 3 ]]$e1c, data = smallData ) )
print( predict( fit3sls[[ 3 ]]$e1c$eq[[ 1 ]], data = smallData ) )

print( predict( fit3sls[[ 4 ]]$e2e, se.fit = TRUE, level = 0.9,
   data = smallData ) )
print( predict( fit3sls[[ 5 ]]$e2e$eq[[ 1 ]], se.pred = TRUE, level = 0.99,
   data = smallData ) )

print( predict( fit3sls[[ 1]]$e3, interval = "prediction", level = 0.975,
   data = smallData ) )
print( predict( fit3sls[[ 1 ]]$e3$eq[[ 1 ]], interval = "confidence", level = 0.8,
   data = smallData ) )

print( predict( fit3sls[[ 2 ]]$e4e, se.fit = TRUE, interval = "confidence",
   level = 0.999, data = smallData ) )
print( predict( fit3sls[[ 2 ]]$e4e$eq[[ 2 ]], se.pred = TRUE, interval = "prediction",
   level = 0.75, data = smallData ) )

print( predict( fit3sls[[ 3 ]]$e5, se.fit = TRUE, interval = "prediction",
   data = smallData ) )
print( predict( fit3sls[[ 3 ]]$e5$eq[[ 1 ]], se.pred = TRUE, interval = "confidence",
   data = smallData ) )

print( predict( fit3slsi[[ 4 ]]$e3e, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, data = smallData ) )
print( predict( fit3slsd[[ 5 ]]$e4$eq[[ 1 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.25, data = smallData ) )


## ************ correlation of predicted values ***************
print( correlation.systemfit( fit3sls[[ 1 ]]$e1c, 2, 1 ) )

print( correlation.systemfit( fit3sls[[ 2 ]]$e2e, 1, 2 ) )

print( correlation.systemfit( fit3sls[[ 3 ]]$e3, 2, 1 ) )

print( correlation.systemfit( fit3sls[[ 4 ]]$e4e, 1, 2 ) )

print( correlation.systemfit( fit3sls[[ 5 ]]$e5, 2, 1 ) )

print( correlation.systemfit( fit3slsi[[ 2 ]]$e3e, 1, 2 ) )

print( correlation.systemfit( fit3slsd[[ 3 ]]$e4, 2, 1 ) )


## ************** F tests ****************
# testing first restriction
print( ftest.systemfit( fit3sls[[ 1 ]]$e1, restrm ) )
print( ftest.systemfit( fit3sls[[ 2 ]]$e1e, restrm ) )
print( ftest.systemfit( fit3sls[[ 3 ]]$e1c, restrm ) )
print( ftest.systemfit( fit3slsi[[ 4 ]]$e1, restrm ) )
print( ftest.systemfit( fit3slsd[[ 5 ]]$e1e, restrm ) )

# testing second restriction
restrOnly2m <- matrix(0,1,7)
restrOnly2q <- 0.5
restrOnly2m[1,2] <- -1
restrOnly2m[1,5] <-  1
# first restriction not imposed 
print( ftest.systemfit( fit3sls[[ 5 ]]$e1c, restrOnly2m, restrOnly2q ) )
print( ftest.systemfit( fit3slsi[[ 1 ]]$e1e, restrOnly2m, restrOnly2q ) )
print( ftest.systemfit( fit3slsd[[ 2 ]]$e1, restrOnly2m, restrOnly2q ) )
# first restriction imposed
print( ftest.systemfit( fit3sls[[ 4 ]]$e2, restrOnly2m, restrOnly2q ) )
print( ftest.systemfit( fit3sls[[ 4 ]]$e3, restrOnly2m, restrOnly2q ) )
print( ftest.systemfit( fit3slsi[[ 5 ]]$e2e, restrOnly2m, restrOnly2q ) )
print( ftest.systemfit( fit3slsi[[ 5 ]]$e3e, restrOnly2m, restrOnly2q ) )
print( ftest.systemfit( fit3slsd[[ 1 ]]$e2, restrOnly2m, restrOnly2q ) )
print( ftest.systemfit( fit3slsd[[ 1 ]]$e3, restrOnly2m, restrOnly2q ) )

# testing both of the restrictions
print( ftest.systemfit( fit3sls[[ 2 ]]$e1e, restr2m, restr2q ) )
print( ftest.systemfit( fit3slsi[[ 3 ]]$e1, restr2m, restr2q ) )
print( ftest.systemfit( fit3slsd[[ 4 ]]$e1e, restr2m, restr2q ) )


## ************** Wald tests ****************
# testing first restriction
print( waldtest.systemfit( fit3sls[[ 1 ]]$e1, restrm ) )
print( waldtest.systemfit( fit3sls[[ 2 ]]$e1e, restrm ) )
print( waldtest.systemfit( fit3sls[[ 3 ]]$e1c, restrm ) )
print( waldtest.systemfit( fit3slsi[[ 4 ]]$e1, restrm ) )
print( waldtest.systemfit( fit3slsd[[ 5 ]]$e1e, restrm ) )

# testing second restriction
# first restriction not imposed
print( waldtest.systemfit( fit3sls[[ 5 ]]$e1c, restrOnly2m, restrOnly2q ) )
print( waldtest.systemfit( fit3slsi[[ 1 ]]$e1e, restrOnly2m, restrOnly2q ) )
print( waldtest.systemfit( fit3slsd[[ 2 ]]$e1, restrOnly2m, restrOnly2q ) )
# first restriction imposed
print( waldtest.systemfit( fit3sls[[ 4 ]]$e2, restrOnly2m, restrOnly2q ) )
print( waldtest.systemfit( fit3sls[[ 4 ]]$e3, restrOnly2m, restrOnly2q ) )
print( waldtest.systemfit( fit3slsi[[ 5 ]]$e2e, restrOnly2m, restrOnly2q ) )
print( waldtest.systemfit( fit3slsi[[ 5 ]]$e3e, restrOnly2m, restrOnly2q ) )
print( waldtest.systemfit( fit3slsd[[ 1 ]]$e2, restrOnly2m, restrOnly2q ) )
print( waldtest.systemfit( fit3slsd[[ 1 ]]$e3, restrOnly2m, restrOnly2q ) )

# testing both of the restrictions
print( waldtest.systemfit( fit3sls[[ 2 ]]$e1e, restr2m, restr2q ) )
print( waldtest.systemfit( fit3slsi[[ 3 ]]$e1, restr2m, restr2q ) )
print( waldtest.systemfit( fit3slsd[[ 4 ]]$e1e, restr2m, restr2q ) )
