### two stage least squares....
# 

#	$Id: systemfit.R,v 1.1 2002/09/18 08:32:44 hamannj Exp $	


# performs two-stage least squares on the system of equations 
ols.systemfit <- function(	
	      			eqns, 
				instruments,
				eqnlabels,
				data )	
{
	results <- list()
	resulti <- list()

	for(i in 1:length( eqns ) ) 
	{

		# perform the two stage least squares regression
		y	  <- eval( attr( terms( eqns[[i]] ), "variables" )[[2]] )
		x	  <- model.matrix( eqns[[i]] )

		v <- diag( dim( model.matrix( eqns[[i]] )[2] ) )
		# two stage least squares results...
		b	  <- solve( t(x) %*% x ) %*% t(x) %*% y
		resids	  <- y - x %*% b
		n	  <- length( y )
		p	  <- ncol( x )
		s	  <- sum(resids^2)/(n - p)
		dfe	  <- n - p
		se	  <- sqrt( diag( solve( t(x) %*% x ) ) * s )
		t	  <- b/se
		prob	  <- 2.0*(1.0 - pt(abs(t), dfe))
		mse	  <- ( s * dfe / dfe )
		rmse	  <- sqrt( mse ) 
		r2	  <- 1.0 - ((t(resids)%*%resids)/(t(y)%*%y-n*mean(y)^2))
		adjr2	  <- 1.0 - ((n-1)/(n-p))*(1.0-r2)
		covb	  <- solve( t(x) %*% x )

		# build the "return" structure for the 2sls part
		resulti$method	     <- "ols"
		resulti$eqnlabel     <- eqnlabels[[i]]			
		resulti$formula	     <- eqns[[i]]
		resulti$dfe	     <- dfe
		resulti$dfm	     <- n - dfe
		resulti$model.matrix <- model.matrix(eqns[[i]] )
		resulti$model.frame  <- model.frame(eqns[[i]] )
		resulti$instruments  <- inst
		resulti$response     <- y
		resulti$predicted    <- x %*% b
		resulti$residuals    <- resids 
		##resulti$ztzinv	     <- ztzinv
		resulti$v	     <- v
		resulti$b	     <- b
		names(resulti$b)     <- colnames( model.matrix( eqns[[i]] ) )
		resulti$n	     <- n
		resulti$s	     <- s
		resulti$sse	     <- s * dfe
		resulti$mse	     <- mse 
		resulti$rmse	     <- rmse 
		resulti$se	     <- se
		resulti$t	     <- t
		resulti$p	     <- prob
		resulti$r2	     <- r2
		resulti$adjr2	     <- adjr2
		resulti$covb	     <- covb
		class(resulti)	     <- "systemfit.ols"
                results[[i]]	     <- resulti

                
	} 

	class(results)	<- "systemfit.system"
	ols <- results

}


## this function produces a table for a single equation 
## in a system of equations
summary.systemfit.ols <- function(object,...)
{
  summary.systemfit.ols <- object
  summary.systemfit.ols
}

## now print the object that comes from the fits...
print.systemfit.ols <- function( x, digits=6, ... )
{

  object <- x

  save.digits <- unlist(options(digits=digits))
  on.exit(options(digits=save.digits))

    cat("\n")
    cat( paste( attr( object, "class" ), "estimates for", object$eqnlabel, "\n" ) )
    
    cat("Model Formula: ")
    print(object$formula)
    cat("Instruments: ")
    print(object$instruments)
    cat("\n")
    
    Signif <- symnum(object$p, corr = FALSE, na = FALSE,
                     cutpoints = c(0,  .001,.01,.05, .1, 1),
                     symbols   =  c("***","**","*","."," "))
    
    table <- cbind(round( object$b, digits ), 
                   round( object$se, digits ), 
                   round( object$t, digits ), 
                   round( object$p, digits ),
                   Signif)
    
    rownames(table) <- names(object$b)
    colnames(table) <- c("Estimate","Std. Error","t value","Pr(>|t|)","")
    
    print.matrix(table, quote = FALSE, right = TRUE )
    cat("---\nSignif. codes: ",attr(Signif,"legend"),"\n")
    
    cat(paste("\nResidual standard error:", round(object$s, digits),
              "on", object$dfe, "degrees of freedom\n"))
    
    cat( paste( "DF-Error:", round(object$dfe, digits),
               "DF-Model:", round(object$dfm, digits),
               "\n" ) )
    
    cat( paste( "SSE:", round(object$sse, digits),
               "MSE:", round(object$s, digits),
               "Root MSE:", round( sqrt(object$s), digits), "\n" ) )
    
    cat( paste( "Multiple R-Squared:", round(object$r2, digits), 
               "Adjusted R-Squared:", round(object$adjr2, digits),
               "\n" ) )
    cat("\n")
  
}



# performs two-stage least squares on the system of equations 
twostage.systemfit <- function(	
			eqns, 
			instruments,
			eqnlabels,
			data )	
{

	results <- list()
	resulti <- list()

	for(i in 1:length( eqns ) ) 
	{

		# perform the two stage least squares regression
		y	  <- eval( attr( terms( eqns[[i]] ), "variables" )[[2]] )
		x	  <- model.matrix( eqns[[i]] )
		z	  <- model.matrix( inst )
		ztzinv	  <- solve( t(z) %*% z )
		v	  <- solve( t(x) %*% z %*% ztzinv %*% t(z) %*% x )

		# two stage least squares results...
		b	  <- v %*% t(x) %*% z %*% ztzinv %*% t(z) %*% y
		resids	  <- y - x %*% b
		n	  <- length( y )
		p	  <- ncol( x )
		s	  <- sum(resids^2)/(n - p)
		dfe	  <- n - p
		se	  <- sqrt(diag(s*v))
		t	  <- b/se
		prob	  <- 2.0*(1.0 - pt(abs(t), dfe))
		mse	  <- ( s * dfe / dfe )
		rmse	  <- sqrt( mse ) 
		r2	  <- 1.0 - ((t(resids)%*%resids)/(t(y)%*%y-n*mean(y)^2))
		adjr2	  <- 1.0 - ((n-1)/(n-p))*(1.0-r2)
		covb	  <- v * s

		# get the residuals from the 2sls on the instruments
		instres <- lsfit( model.frame( inst ), model.matrix( eqns[[i]] ) )$coef
		temp2 <- model.matrix( inst ) %*% instres
		resulti$instres <- instres
		resulti$tslsres <- temp2

		# build the "return" structure for the 2sls part
		resulti$method	     <- "2sls"
		resulti$eqnlabel     <- eqnlabels[[i]]			
		resulti$formula	     <- eqns[[i]]
		resulti$dfe	     <- dfe
		resulti$dfm	     <- n - dfe
		resulti$model.matrix <- model.matrix(eqns[[i]] )
		resulti$model.frame  <- model.frame(eqns[[i]] )
		resulti$instruments  <- inst
		resulti$response     <- y
		resulti$predicted    <- x %*% b
		resulti$residuals    <- resids 
		resulti$ztzinv	     <- ztzinv
		resulti$v	     <- v
		resulti$b	     <- b
		names(resulti$b)     <- colnames( model.matrix( eqns[[i]] ) )
		resulti$n	     <- n
		resulti$s	     <- s
		resulti$sse	     <- s * dfe
		resulti$mse	     <- mse 
		resulti$rmse	     <- rmse 
		resulti$se	     <- se
		resulti$t	     <- t
		resulti$p	     <-  prob
		resulti$r2	     <- r2
		resulti$adjr2	     <- adjr2
		resulti$covb	     <- covb
		class(resulti)	     <- "systemfit.twostage"
		results[[i]]	     <- resulti

	} 

	class(results)	<- "systemfit.system"
	twostage <- results
}


summary.systemfit.twostage <- function(object,...)
{
  summary.systemfit.twostage <- object
  summary.systemfit.twostage
}

## now print the object that comes from the fits...
print.systemfit.twostage <- function( x,digits=6,... )
{

  object <- x

  save.digits <- unlist(options(digits=digits))
  on.exit(options(digits=save.digits))

    cat("\n")
    cat( paste( attr( object, "class" ), "estimates for", object$eqnlabel, "\n" ) )
    
    cat("Model Formula: ")
    print(object$formula)
    cat("Instruments: ")
    print(object$instruments)
    cat("\n")
    
    Signif <- symnum(object$p, corr = FALSE, na = FALSE,
                     cutpoints = c(0,  .001,.01,.05, .1, 1),
                     symbols   =  c("***","**","*","."," "))
    
    table <- cbind(round( object$b, digits ), 
                   round( object$se, digits ), 
                   round( object$t, digits ), 
                   round( object$p, digits ),
                   Signif)
    
    rownames(table) <- names(object$b)
    colnames(table) <- c("Estimate","Std. Error","t value","Pr(>|t|)","")
    
    print.matrix(table, quote = FALSE, right = TRUE )
    cat("---\nSignif. codes: ",attr(Signif,"legend"),"\n")
    
    cat(paste("\nResidual standard error:", round(object$s, digits),
              "on", object$dfe, "degrees of freedom\n"))
    
    cat( paste( "DF-Error:", round(object$dfe, digits),
               "DF-Model:", round(object$dfm, digits),
               "\n" ) )
    
    cat( paste( "SSE:", round(object$sse, digits),
               "MSE:", round(object$s, digits),
               "Root MSE:", round( sqrt(object$s), digits), "\n" ) )
    
    cat( paste( "Multiple R-Squared:", round(object$r2, digits), 
               "Adjusted R-Squared:", round(object$adjr2, digits),
               "\n" ) )
    cat("\n")
  
}


# add the SUR firs here...
sur.systemfit <- function(	
                          eqns, 
                          instruments,
                          eqnlabels,
                          data )	
{

	results <- list()
	resulti <- list()
	u2	<-  matrix( 0, dim(data)[1], length( eqns ) )

	# perform the two-stage fits
	osls <- ols.systemfit(	
                              eqns, 
                              instruments,
                              eqnlabels,
                              data )	
        

	# these are the ones that wil be used to build the big matrix
	t3 <- NULL
	bigb <- NULL
	bigy <- NULL
	bigt <- NULL
	bigse <- NULL
	bigp <- NULL

	for(i in 1:length( eqns ) ) 
	{
		# build the final large matrix...
		tr <- NULL
	
		# get the dimensions of the current matrix
		for(j in 1:length( eqns ) ) 
		{
			if( i == j )
			{
				##tr <- cbind( tr, tsls[[i]]$tslsres )
				##tr <- cbind( tr, osls[[i]]$residuals )
				tr <- cbind( tr, osls[[i]]$model.matrix )
			}		
			else
			{
				# bind the zero matrix to the row
				di <- dim( model.matrix( eqns[[j]] ) )[1]
				dj <- dim( model.matrix( eqns[[j]] ) )[2]
				tr <- cbind( tr, matrix( 0, di, dj ) )
			}
		}

		t3 <- rbind( t3, tr )

		# now add the rows to the bigX matrix
		# or should this be the new fitted y values
		# from the two stage least squares fits...
		y      <- eval( attr( terms( eqns[[i]] ), "variables" )[[2]] )
		bigy <- rbind( bigy, matrix( y ) )
	} 

	## get the variance-covariance matrix from the two stage results
	varcov <- varcov.systemfit( osls )	

	parta <- kronecker( solve( varcov ), 
			    diag( dim( model.matrix( eqns[[1]] ) )[1] ) )
	part1 <- solve( t(t3) %*% parta %*% t3 ) # covariance matrix
	part2 <- t(t3) %*% parta %*% bigy
	bigb <- part1 %*% part2

	# compute the se, t, and p values...
	bigse <- matrix( sqrt( diag( part1 ) ) )
	bigt <- bigb/bigse

	# extract the results
	idx <- matrix( 0, length( eqns ), 2 )
	for(i in 1:length( eqns ) ) 
	{


		## get the index for stripping out the estimates
		if( i == 1 )
		{
			idx[i,1] <- 1
			idx[i,2] <- dim( model.matrix( eqns[[i]] ) )[2]
		}      
		else
		{
			idx[i,1] <- idx[i-1,2]+1
			idx[i,2] <- idx[i,1] + 
				    dim( model.matrix(eqns[[i]]))[2]-1	
		}

		start1	<- idx[i,1]
		start2	<- idx[i,2]		

		# tree stage least squares results...
		x	  <- model.matrix( eqns[[i]] )
		y	  <- eval( attr( terms( eqns[[i]] ), "variables" )[[2]] )
		b	  <- matrix( bigb[start1:start2] )
		resids	  <- y - x %*% b
		n	  <- length( y )
		p	  <- ncol( x )
		s	  <- sum(resids^2)/(n - p)
		dfe	  <- n - p
		se	  <- matrix( bigse[start1:start2] )
		t	  <- matrix( bigt[start1:start2] )
		prob	  <- 2.0*(1.0 - pt(abs(t), dfe))
		mse	  <- ( s * dfe / dfe )
		rmse	  <- sqrt( mse ) 
		r2	  <- 1.0 - ((t(resids)%*%resids)/(t(y)%*%y-n*mean(y)^2))
		adjr2	  <- 1.0 - ((n-1)/(n-p))*(1.0-r2)

		# get the parameter var-cov matrix fo the eq
		icol	  <- ncol( model.matrix( eqns[[i]] ) )
		jcol	  <- ncol( model.matrix( eqns[[i]] ) )
		startrow  <- idx[i,1]
		endrow	  <- idx[i,2]
		startcol  <- idx[i,1]
		endcol	  <- idx[i,2]
		covb	  <- matrix( part1[startrow:endrow,startcol:endcol], icol, jcol )

		# build the "return" structure for the 3sls part
		resulti$method	     <- "sur"
		resulti$eqnlabel     <- eqnlabels[[i]]			
		resulti$formula	     <- eqns[[i]]
		resulti$dfe	     <- dfe
		resulti$dfm	     <- n - dfe
		resulti$model.matrix <- model.matrix(eqns[[i]] )
		resulti$model.frame  <- model.frame(eqns[[i]] )
		resulti$instruments  <- inst
		resulti$response     <- osls[[i]]$reponse
		resulti$predicted    <- x %*% b
		resulti$residuals    <- resids 
		##resulti$ztzinv	     <- tsls[[i]]$ztzinv
		resulti$v	     <- osls[[i]]$v
		resulti$b	     <- b
		names(resulti$b)     <- colnames( model.matrix( eqns[[i]] ) )
		resulti$n	     <- n
		resulti$s	     <- s
		resulti$sse	     <- s * dfe
		resulti$mse	     <- mse 
		resulti$rmse	     <- rmse 
		resulti$se	     <- se
		resulti$t	     <- t
		resulti$p	     <- prob
		resulti$r2	     <- r2
		resulti$adjr2	     <- adjr2
		resulti$systemcovb   <- part1	# sneak in the whole covariance matrix 
		resulti$covb	     <- covb	# this is the var-cov matrix for the eq  
		class(resulti)	     <- "systemfit.sur"
		results[[i]]	     <- resulti

	}

	class(results)	<- "systemfit.system"
	sur <- results
}


summary.systemfit.sur <- function(object,...)
{
  summary.systemfit.sur <- object
  summary.systemfit.sur
}

## now print the object that comes from the fits...
print.systemfit.sur <- function( x,digits=6,... )
{

  object <- x

  save.digits <- unlist(options(digits=digits))
  on.exit(options(digits=save.digits))

    cat("\n")
    cat( paste( attr( object, "class" ), "estimates for", object$eqnlabel, "\n" ) )
    
    cat("Model Formula: ")
    print(object$formula)
    cat("Instruments: ")
    print(object$instruments)
    cat("\n")
    
    Signif <- symnum(object$p, corr = FALSE, na = FALSE,
                     cutpoints = c(0,  .001,.01,.05, .1, 1),
                     symbols   =  c("***","**","*","."," "))
    
    table <- cbind(round( object$b, digits ), 
                   round( object$se, digits ), 
                   round( object$t, digits ), 
                   round( object$p, digits ),
                   Signif)
    
    rownames(table) <- names(object$b)
    colnames(table) <- c("Estimate","Std. Error","t value","Pr(>|t|)","")
    
    print.matrix(table, quote = FALSE, right = TRUE )
    cat("---\nSignif. codes: ",attr(Signif,"legend"),"\n")
    
    cat(paste("\nResidual standard error:", round(object$s, digits),
              "on", object$dfe, "degrees of freedom\n"))
    
    cat( paste( "DF-Error:", round(object$dfe, digits),
               "DF-Model:", round(object$dfm, digits),
               "\n" ) )
    
    cat( paste( "SSE:", round(object$sse, digits),
               "MSE:", round(object$s, digits),
               "Root MSE:", round( sqrt(object$s), digits), "\n" ) )
    
    cat( paste( "Multiple R-Squared:", round(object$r2, digits), 
               "Adjusted R-Squared:", round(object$adjr2, digits),
               "\n" ) )
    cat("\n")
  
}





# performs three-stage least squares on the system of equations 
threestage.systemfit <- function(	
					eqns, 
					instruments,
					eqnlabels,
					data )	
{

	results <- list()
	resulti <- list()
	u2	<-  matrix( 0, dim(data)[1], length( eqns ) )

	# perform the two-stage fits
	tsls <- twostage.systemfit(	
					eqns, 
					instruments,
					eqnlabels,
					data )	


	# these are the ones that wil be used to build the big matrix
	t3 <- NULL
	bigb <- NULL
	bigy <- NULL
	bigt <- NULL
	bigse <- NULL
	bigp <- NULL

	for(i in 1:length( eqns ) ) 
	{
		# build the final large matrix...
		tr <- NULL
	
		# get the dimensions of the current matrix
		for(j in 1:length( eqns ) ) 
		{
			if( i == j )
			{
				tr <- cbind( tr, tsls[[i]]$tslsres )
			}		
			else
			{
				# bind the zero matrix to the row
				di <- dim( model.matrix( eqns[[j]] ) )[1]
				dj <- dim( model.matrix( eqns[[j]] ) )[2]
				tr <- cbind( tr, matrix( 0, di, dj ) )
			}
		}

		t3 <- rbind( t3, tr )

		# now add the rows to the bigX matrix
		# or should this be the new fitted y values
		# from the two stage least squares fits...
		y      <- eval( attr( terms( eqns[[i]] ), "variables" )[[2]] )
		bigy <- rbind( bigy, matrix( y ) )
	} 

	## get the variance-covariance matrix from the two stage results
	varcov <- varcov.systemfit( tsls )	

	parta <- kronecker( solve( varcov ), 
			    diag( dim( model.matrix( eqns[[1]] ) )[1] ) )
	part1 <- solve( t(t3) %*% parta %*% t3 ) # covariance matrix
	part2 <- t(t3) %*% parta %*% bigy
	bigb <- part1 %*% part2

	# compute the se, t, and p values...
	bigse <- matrix( sqrt( diag( part1 ) ) )
	bigt <- bigb/bigse

	# extract the results
	idx <- matrix( 0, length( eqns ), 2 )
	for(i in 1:length( eqns ) ) 
	{


		## get the index for stripping out the estimates
		if( i == 1 )
		{
			idx[i,1] <- 1
			idx[i,2] <- dim( model.matrix( eqns[[i]] ) )[2]
		}      
		else
		{
			idx[i,1] <- idx[i-1,2]+1
			idx[i,2] <- idx[i,1] + 
				    dim( model.matrix(eqns[[i]]))[2]-1	
		}

		start1	<- idx[i,1]
		start2	<- idx[i,2]		

		# tree stage least squares results...
		x	  <- model.matrix( eqns[[i]] )
		y	  <- eval( attr( terms( eqns[[i]] ), "variables" )[[2]] )
		b	  <- matrix( bigb[start1:start2] )
		resids	  <- y - x %*% b
		n	  <- length( y )
		p	  <- ncol( x )
		s	  <- sum(resids^2)/(n - p)
		dfe	  <- n - p
		se	  <- matrix( bigse[start1:start2] )
		t	  <- matrix( bigt[start1:start2] )
		prob	  <- 2.0*(1.0 - pt(abs(t), dfe))
		mse	  <- ( s * dfe / dfe )
		rmse	  <- sqrt( mse ) 
		r2	  <- 1.0 - ((t(resids)%*%resids)/(t(y)%*%y-n*mean(y)^2))
		adjr2	  <- 1.0 - ((n-1)/(n-p))*(1.0-r2)

		# get the parameter var-cov matrix fo the eq
		icol	  <- ncol( model.matrix( eqns[[i]] ) )
		jcol	  <- ncol( model.matrix( eqns[[i]] ) )
		startrow  <- idx[i,1]
		endrow	  <- idx[i,2]
		startcol  <- idx[i,1]
		endcol	  <- idx[i,2]
		covb	  <- matrix( part1[startrow:endrow,startcol:endcol], icol, jcol )

		# build the "return" structure for the 3sls part
		resulti$method	     <- "3sls"
		resulti$eqnlabel     <- eqnlabels[[i]]			
		resulti$formula	     <- eqns[[i]]
		resulti$dfe	     <- dfe
		resulti$dfm	     <- n - dfe
		resulti$model.matrix <- model.matrix(eqns[[i]] )
		resulti$model.frame  <- model.frame(eqns[[i]] )
		resulti$instruments  <- inst
		resulti$response     <- tsls[[i]]$reponse
		resulti$predicted    <- x %*% b
		resulti$residuals    <- resids 
		resulti$ztzinv	     <- tsls[[i]]$ztzinv
		resulti$v	     <- tsls[[i]]$v
		resulti$b	     <- b
		names(resulti$b)     <- colnames( model.matrix( eqns[[i]] ) )
		resulti$n	     <- n
		resulti$s	     <- s
		resulti$sse	     <- s * dfe
		resulti$mse	     <- mse 
		resulti$rmse	     <- rmse 
		resulti$se	     <- se
		resulti$t	     <- t
		resulti$p	     <- prob
		resulti$r2	     <- r2
		resulti$adjr2	     <- adjr2
		resulti$systemcovb   <- part1	# sneak in the whole covariance matrix 
		resulti$covb	     <- covb	# this is the var-cov matrix for the eq  
		class(resulti)	     <- "systemfit.threestage"
		results[[i]]	     <- resulti

	}

	class(results)	<- "systemfit.system"
	threestage <- results
}



summary.systemfit.threestage <- function(object,...)
{
  summary.systemfit.threestage <- object
  summary.systemfit.threestage
}

## now print the object that comes from the fits...
print.systemfit.threestage <- function( x,digits=6,... )
{

  object <- x
  
  save.digits <- unlist(options(digits=digits))
  on.exit(options(digits=save.digits))

    cat("\n")
    cat( paste( attr( object, "class" ), "estimates for", object$eqnlabel, "\n" ) )
    
    cat("Model Formula: ")
    print(object$formula)
    cat("Instruments: ")
    print(object$instruments)
    cat("\n")
    
    Signif <- symnum(object$p, corr = FALSE, na = FALSE,
                     cutpoints = c(0,  .001,.01,.05, .1, 1),
                     symbols   =  c("***","**","*","."," "))
    
    table <- cbind(round( object$b, digits ), 
                   round( object$se, digits ), 
                   round( object$t, digits ), 
                   round( object$p, digits ),
                   Signif)
    
    rownames(table) <- names(object$b)
    colnames(table) <- c("Estimate","Std. Error","t value","Pr(>|t|)","")
    
    print.matrix(table, quote = FALSE, right = TRUE )
    cat("---\nSignif. codes: ",attr(Signif,"legend"),"\n")
    
    cat(paste("\nResidual standard error:", round(object$s, digits),
              "on", object$dfe, "degrees of freedom\n"))
    
    cat( paste( "DF-Error:", round(object$dfe, digits),
               "DF-Model:", round(object$dfm, digits),
               "\n" ) )
    
    cat( paste( "SSE:", round(object$sse, digits),
               "MSE:", round(object$s, digits),
               "Root MSE:", round( sqrt(object$s), digits), "\n" ) )
    
    cat( paste( "Multiple R-Squared:", round(object$r2, digits), 
               "Adjusted R-Squared:", round(object$adjr2, digits),
               "\n" ) )
    cat("\n")
  
}


## this function returns the variance-covariance matrix
## from the results set for equation ij
threestage.cov <- function( results, eqni, eqnj )
{

	## get the information about eqni and enqj
	## get the size of the array for the matrix 
	## you are going to extract 
	icol	  <- ncol( results[[eqni]]$model.matrix )
	jcol	  <- ncol( results[[eqnj]]$model.matrix )
	
	## now get the offsets
	## 1 - start row
	## 2 - end row

	rows <- matrix( 0, length( results ), 2 )
	#cols <- matrix( 0, length( results ), 2 )
	for(i in 1:length( results ) ) 
	{
		## get the index for stripping out the estimates
		if( i == 1 )
		{
			rows[i,1] <- 1
			rows[i,2] <- dim( results[[i]]$model.matrix)[2]
		}      
		else
		{
			rows[i,1] <- rows[i-1,2]+1
			rows[i,2] <- rows[i,1] + dim( results[[i]]$model.matrix )[2]-1
		}

	}		

	startrow <- rows[eqni,1]
	endrow	 <- rows[eqni,2]
	
	startcol <- rows[eqnj,1]
	endcol	 <- rows[eqnj,2]

	test <- matrix( results[[1]]$systemcovb[startrow:endrow,startcol:endcol], icol, jcol )
		
}

## this function returns test statistic for
## the hausman test which.... i forget, but people want to see it... 
## from the sas docs
## given 2 estimators, b0 abd b1, where under the null hypothesis,
## both are consistent, but only b0 is asympt. efficient and 
## under the alter. hypo only b1 is consistent, so the statistic (m) is
hausman.systemfit <- function( results0, results1 )
{

	v0 <- results0$covb
	v1 <- results1$covb
	q  <- results1$b - results0$b
		
	hausman <- t( q ) %*% ( v1 - v0 ) %*% q

}


## this function returns the covariance of the residuals
## the method will return the same matrix values as are 
## returned in SAS in proc model
varcov.systemfit <- function( results )	
{

	u2 <-  matrix( 0, length( results ), length( results ) )

	## use bind to create a vector for the residuals 
	for(i in 1:length( results ) ) 
	{	
		## use bind to create a vector for the residuals 
		for(j in 1:length( results ) ) 
		{	
			ri	<- results[[i]]$residuals
			dfei	<- results[[i]]$dfe

			rj	<- results[[j]]$residuals
			dfej	<- results[[j]]$dfe

			# from SAS
			cvij <- ( t( ri ) %*% rj ) / ( sqrt( dfei * dfej ) )
		
			u2[i,j] <- cvij
		}

	}

	varcov <- u2
	varcov
}





  
  ## this function returns a vector of the 
  ## cross-equation corrlations between eq i and eq j
  ## from the results set for equation ij
correlation.systemfit <- function( results, eqni, eqnj )
{

	cij <- threestage.cov( results, eqni, eqnj )
	cii <- threestage.cov( results, eqni, eqni )
	cjj <- threestage.cov( results, eqnj, eqnj )

	rij <- NULL

	for(i in 1:results[[1]]$n ) 
	{

		xik <- model.matrix( results[[eqni]]$formula )[i,]
		xjk <- model.matrix( results[[eqnj]]$formula )[i,]

		top <- xik %*% cij %*% xjk
		bottom <- sqrt( ( xik %*% cii %*% xik ) * ( xjk %*% cjj %*% xjk ) )
		rijk <- top / bottom

		rij <- rbind( rij, rijk )
	
	}

	correlation <- rij  
	correlation
}


## this function returns a vector of the 
## cross-equation corrlations between eq i and eq j
## from the results set for equation ij
## you need to put some check in here to make sure both
## are the name type

## determines the improvement of resultsj (3sls) over 
## resultsi (2sls) for equation i and returns a matrix
## of the values, so you can examine the range, mean, etc
se.ratio.systemfit <- function( resultsi, resultsj, eqni )
{

	ratio <- NULL

	for(i in 1:resultsi[[1]]$n ) 
	{

		xik <- model.matrix( resultsi[[eqni]]$formula )[i,]

		top    <- sqrt( xik %*% resultsi[[eqni]]$covb %*% xik )
		bottom <- sqrt( xik %*% resultsj[[eqni]]$covb %*% xik )
		rk     <- top / bottom

		ratio <- rbind( ratio, rk )
	
	}

	se.ratio <- ratio  
	se.ratio
}

summary.systemfit.system <- function(object,...)
{
  summary.systemfit.system <- object
  summary.systemfit.system
}


## now print the object that comes from the fits...
print.systemfit.system <- function( x,digits=6,... )
{

  object <- x
  
  save.digits <- unlist(options(digits=digits))
  on.exit(options(digits=save.digits))


  table <- NULL
  labels <- NULL
  
  cat("\n")
  cat("systemfit results: \n")
  cat("\n")

  for(i in 1:length( object ) ) 
    {
      row <- NULL
      row <- cbind(
                   round( object[[i]]$dfm, digits ),
                   round( object[[i]]$dfe, digits ), 
                   round( object[[i]]$sse, digits ),
                   round( object[[i]]$mse, digits ),
                   round( object[[i]]$rmse, digits ),
                   round( object[[i]]$r2, digits ),
                   round( object[[i]]$adjr2, digits ) )
      
      table <- rbind( table, row )
      labels <- rbind( labels, object[[i]]$eqnlabel ) 
    }
  
  rownames(table) <- c( labels )
  colnames(table) <- c(	
                       "DF Model",
                       "DF Error",
                       "SSE",
                       "MSE",
                       "RMSE",
                       "R2",
                       "Adj R2" )
  
  print.matrix(table, quote = FALSE, right = TRUE )
  cat("\n")

  ### if the system is sur/3sls then print the
  ### varcov matrix used for estimation
  ###cat( paste( attr( object, "class" ), "estimates for", object$eqnlabel, "\n" ) )
  
  cat("The variance-covariance matrix\n")
  vc <- varcov.systemfit( object )
  rownames(vc) <- labels
  colnames(vc) <- labels
  print( vc )

  
  ## now print the individual equations
  for(i in 1:length( object ) ) 
    {
      print( object[[i]], digits )
    }

    save.digits <- unlist(options(digits=digits))
    on.exit(options(digits=save.digits))
  
}
