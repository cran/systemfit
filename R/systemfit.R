### this is the systemfit project....

#	$Id: systemfit.R,v 1.18 2003/11/10 20:19:04 hamannj Exp $	

systemfit <- function( method,
                        eqns,
                        eqnlabels=c(as.character(1:length(eqns))),
                        inst=NULL,
                        data=list(),
                        R.restr=NULL,
                        q.restr=matrix(0,max(nrow(R.restr),0),1),
                        TX=NULL,
                        maxiter=1,
                        tol=1e-5,
                        rcovformula=1,
                        formula3sls="GLS",
                        probdfsys=!(is.null(R.restr) & is.null(TX)),
                        single.eq.sigma=(is.null(R.restr) & is.null(TX)),
                        solvetol=.Machine$double.eps)
{
  attach(data)

   ## some tests
   if(!(method=="OLS" | method=="WLS" | method=="SUR" | method=="2SLS" | method=="W2SLS" |
      method=="3SLS")){
      stop("The method must be 'OLS', 'WLS', 'SUR', '2SLS', 'W2SLS' or '3SLS'")}
   if((method=="2SLS" | method=="W2SLS" | method=="3SLS") & is.null(inst)) {
      stop("The methods '2SLS', 'W2SLS' and '3SLS' need instruments!")}


  results <- list()               # results to be returned
  results$eq <- list()            # results for the individual equations
  resulti <- list()               # results of the ith equation
  residi  <- list()               # residuals equation wise
  iter    <- NULL                 # number of iterations
  G       <- length( eqns )       # number of equations
  y       <- list()               # endogenous variables equation wise
  Y       <- matrix( 0, 0, 1 )    # stacked endogenous variables
  x       <- list()               # regressors equation-wise
  X       <- matrix( 0, 0, 0 )    # stacked matrices of all regressors (unrestricted)
  n       <- array( 0, c(G))      # number of observations in each equation
  k       <- array( 0, c(G) )     # number of (unrestricted) coefficients/regressors in each equation
  instl   <- list()               # list of the instruments for each equation
  ssr     <- array( 0, c(G))      # sum of squared residuals of each equation
  mse     <- array( 0, c(G))      # mean square error (residuals) of each equation
  rmse    <- array( 0, c(G))      # root of mse
  r2      <- array( 0, c(G))      # R-squared value
  adjr2   <- array( 0, c(G))      # adjusted R-squared value
  xnames  <- NULL                 # names of regressors

  for(i in 1:G )  {
    y[[i]] <-  eval( attr( terms( eqns[[i]] ), "variables" )[[2]] )
    Y      <-  c( Y, y[[i]] )
    x[[i]] <-  model.matrix( eqns[[i]] )
    X      <-  rbind( cbind( X, matrix( 0, nrow( X ), ncol( x[[i]] ))),
                       cbind( matrix( 0, nrow( x[[i]] ), ncol( X )), x[[i]]))
    n[i]   <-  length( y[[i]] )
    k[i]   <-  ncol(x[[i]])
    for(j in 1:k[i]) {
      xnames <- c( xnames, paste("eq",as.character(i),colnames( x[[i]] )[j] ))
    }
  }

  N  <- sum( n )              # total number of observations
  K  <- sum( k )              # total number of (unrestricted) coefficients/regressors
  Ki <- K                     # total number of linear independent coefficients
  ki <- k                     # total number of linear independent coefficients in each equation
  if(!is.null(TX)) {
    XU <- X
    X  <- XU %*% TX
    Ki <- Ki - ( nrow( TX ) - ncol( TX ) )
    for(i in 1:G) {
       ki[i] <- ncol(X)
       for(j in 1: ncol(X) ) {
          if(sum(X[(1+sum(n[1:i])-n[i]):(sum(n[1:i])),j]^2)==0) ki[i] <- ki[i]-1
       }
    }
  }
  if(!is.null(R.restr)) {
    Ki  <- Ki - nrow(R.restr)
    if(is.null(TX)) {
       for(j in 1:nrow(R.restr)) {
          for(i in 1:G) {                   # search for restrictions that are NOT cross-equation
             if(sum(R.restr[j,(1+sum(k[1:i])-k[i]):(sum(k[1:i]))]^2)==sum(R.restr[j,]^2)) {
                ki[i] <- ki[i]-1
             }
          }
       }
    }
  }
  df <- n - ki              # degress of freedom of each equation

  ## only for OLS, WLS and SUR estimation
  if(method=="OLS" | method=="WLS" | method=="SUR") {
    if(is.null(R.restr)) {
      b <- solve( crossprod( X ), crossprod( X, Y ), tol=solvetol )  # estimated coefficients
    } else {
      W <- rbind( cbind( t(X) %*% X, t(R.restr) ),
                  cbind( R.restr, matrix( 0, nrow(R.restr), nrow(R.restr) )))
      V <- rbind( t(X) %*% Y , q.restr )
      b <- ( solve( W, tol=solvetol ) %*% V )[1:ncol(X)]
    }
  }

  ## only for OLS estimation
  if(method=="OLS") {
    resids <- Y - X %*% b                                        # residuals
    for(i in 1:G) residi[[i]] <- resids[(1+sum(n[1:i])-n[i]):(sum(n[1:i]))]
    if(single.eq.sigma) {
      rcov  <- matrix( 0, G, G )                                 # residual covariance matrix
      for(i in 1:G) rcov[i,i] <- sum(residi[[i]]*residi[[i]])/df[i]
      Oinv   <- solve( rcov, tol=solvetol ) %x% diag(1,n[1],n[1])              # Omega inverse
      if(is.null(R.restr)) {
        bcov   <- solve( t(X) %*% Oinv %*% X, tol=solvetol )    # coefficient covariance matrix
      } else {
        W <- rbind( cbind( t(X) %*% Oinv %*% X, t(R.restr) ),
                    cbind( R.restr, matrix( 0, nrow(R.restr), nrow(R.restr) )))
        bcov <- solve( W, tol=solvetol )[1:ncol(X),1:ncol(X)]
      }
    } else {
      s2     <- sum(resids^2)/(N-Ki)                            # sigma squared
      if(is.null(R.restr)) {
        bcov   <- s2 * solve( crossprod( X ), tol=solvetol )    # coefficient covariance matrix
      } else {
        bcov   <- s2 * solve( W, tol=solvetol )[1:ncol(X),1:ncol(X)]  # coefficient covariance matrix
      }
    }
  }

  ## only for WLS estimation
  if(method=="WLS") {
    bl    <- b                       # coefficients of previous step
    bdif  <- b                       # difference of coefficients between this and previous step
    rcov  <- matrix( 0, G, G )                               # residual covariance matrix
    iter  <- 0
    while((sum(bdif^2)/sum(bl^2))^0.5>tol & iter < maxiter) {
      iter  <- iter+1
      bl    <- b                           # coefficients of previous step
      resids <- Y - X %*% b                     # residuals
      for(i in 1:G) residi[[i]] <- resids[(1+sum(n[1:i])-n[i]):(sum(n[1:i]))]
      for(i in 1:G)  {
        if(rcovformula==0) {
          rcov[i,i] <- sum(residi[[i]]*residi[[i]])/n[i]
        } else {
          rcov[i,i] <- sum(residi[[i]]*residi[[i]])/df[i]
        }
      }
      Oinv <- solve( rcov, tol=solvetol ) %x% diag(1,n[1],n[1]) # Omega inverse (= weight. matrix)
      if(is.null(R.restr)) {
        b  <- solve(t(X) %*% Oinv %*% X, tol=solvetol) %*% t(X) %*% Oinv %*%Y   # coefficients
      } else {
        W <- rbind( cbind( t(X) %*% Oinv %*% X, t(R.restr) ),
                    cbind( R.restr, matrix(0, nrow(R.restr), nrow(R.restr) )))
        V <- rbind( t(X) %*% Oinv %*% Y , q.restr )
        Winv <- solve( W, tol=solvetol )
        b <- ( Winv %*% V )[1:ncol(X)]     # restricted coefficients
      }
      bdif <- b-bl                 # difference of coefficients between this and previous step
    }
    if(is.null(R.restr)) {
      bcov <- solve(t(X) %*% Oinv %*% X, tol=solvetol ) # final step coefficient covariance matrix
    } else {
      bcov   <- Winv[1:ncol(X),1:ncol(X)]     # coefficient covariance matrix
    }
    resids <- Y - X %*% b                        # residuals
    for(i in 1:G) residi[[i]] <- resids[(1+sum(n[1:i])-n[i]):(sum(n[1:i]))]
  }

  ## only for SUR estimation
  if(method=="SUR") {
    bl    <- b                       # coefficients of previous step
    bdif  <- b                       # difference of coefficients between this and previous step
    rcov  <- matrix( 0, G, G )                               # residual covariance matrix
    iter  <- 0
    while((sum(bdif^2)/sum(bl^2))^0.5>tol & iter < maxiter) {
      iter  <- iter+1
      bl    <- b                           # coefficients of previous step
      resids <- Y-X%*%b                     # residuals
      for(i in 1:G) residi[[i]] <- resids[(1+sum(n[1:i])-n[i]):(sum(n[1:i]))]
      for(i in 1:G) {
        for(j in 1:G) {
          if(rcovformula==0 | rcovformula==1) {
            rcov[i,j] <- sum(residi[[i]]*residi[[j]])/(
                            sqrt((n[i]-rcovformula*ki[i])*(n[j]-rcovformula*ki[j])))
          } else {
            rcov[i,j] <- sum(residi[[i]]*residi[[j]])/(
                            n[i]-ki[i]-ki[j] + sum( diag(
                            solve( crossprod( x[[i]] ), tol=solvetol ) %*%
                            crossprod( x[[i]], x[[j]]) %*%
                            solve( crossprod( x[[j]] ), tol=solvetol ) %*%
                            crossprod( x[[j]], x[[i]] ))))

          }
        }
      }
      Oinv <- solve( rcov, tol=solvetol ) %x% diag(1,n[1],n[1])  # Omega inverse (= weighting matrix)
      if(is.null(R.restr)) {
        b  <- solve(t(X) %*% Oinv %*% X, tol=solvetol) %*% t(X) %*% Oinv %*%Y   # coefficients
      } else {
        W <- rbind( cbind( t(X) %*% Oinv %*% X, t(R.restr) ),
                    cbind( R.restr, matrix(0, nrow(R.restr), nrow(R.restr) )))
        V <- rbind( t(X) %*% Oinv %*% Y , q.restr )
        Winv <- solve( W, tol=solvetol )
        b <- ( Winv %*% V )[1:ncol(X)]     # restricted coefficients
      }
      bdif <- b-bl                         # difference of coefficients between this and previous step
    }
    if(is.null(R.restr)) {
      bcov <- solve(t(X) %*% Oinv %*% X, tol=solvetol )   # final step coefficient covariance matrix
    } else {
      bcov   <- Winv[1:ncol(X),1:ncol(X)]     # coefficient covariance matrix
    }
    resids <- Y - X %*% b                        # residuals
    for(i in 1:G) residi[[i]] <- resids[(1+sum(n[1:i])-n[i]):(sum(n[1:i]))]
  }

  ## only for 2SLS, W2SLS and 3SLS estimation
  if(method=="2SLS" | method=="W2SLS" | method=="3SLS") {
    for(i in 1:G) {
      if(is.list(inst)) { instl[[i]] <- inst[[i]]
      } else              instl[[i]] <- inst
    }
    Xf <- array(0,c(0,ncol(X)))       # fitted X values
    H  <- matrix( 0, 0, 0 )           # stacked matrices of all instruments
    h  <- list()
    for(i in 1:G) {
      Xi <- X[(1+sum(n[1:i])-n[i]):(sum(n[1:i])),]  # regressors of the ith equation (including zeros)
      h[[i]] <- model.matrix( instl[[i]] )
      Xf <- rbind(Xf, h[[i]] %*% solve( crossprod( h[[i]]) , tol=solvetol )
              %*% crossprod( h[[i]], Xi ))                                # 'fitted' X-values
      H  <-  rbind( cbind( H, matrix( 0, nrow( H ), ncol( h[[i]] ))),
                         cbind( matrix( 0, nrow( h[[i]] ), ncol( H )), h[[i]]))

    }
    if(is.null(R.restr)) {
      b <- solve( crossprod( Xf ), crossprod( Xf, Y ), tol=solvetol )  # 2nd stage coefficients
    } else {
      W <- rbind( cbind( crossprod(Xf), t(R.restr) ),
                  cbind( R.restr, matrix(0, nrow(R.restr), nrow(R.restr))))
      V <- rbind( t(Xf) %*% Y , q.restr )
      b <- ( solve( W, tol=solvetol ) %*% V )[1:ncol(X)]     # restricted coefficients
    }
    b2 <- b
  }

  ## only for 2SLS estimation
  if(method=="2SLS") {
    resids <- Y - X %*% b                        # residuals
    for(i in 1:G) residi[[i]] <- resids[(1+sum(n[1:i])-n[i]):(sum(n[1:i]))]
    if(single.eq.sigma) {
      rcov  <- matrix( 0, G, G )                                 # residual covariance matrix
      for(i in 1:G) rcov[i,i] <- sum(residi[[i]]*residi[[i]])/(df[i])
      Oinv   <- solve( rcov, tol=solvetol ) %x% diag(1,n[1],n[1])              # Omega inverse
      if(is.null(R.restr)) {
        bcov   <- solve( t(Xf) %*% Oinv %*% Xf, tol=solvetol )           # coefficient covariance matrix
      } else {
        W <- rbind( cbind( t(Xf) %*% Oinv %*% Xf, t(R.restr) ),
                    cbind( R.restr, matrix( 0, nrow(R.restr), nrow(R.restr) )))
        bcov <- solve( W, tol=solvetol )[1:ncol(X),1:ncol(X)]
      }
    } else {
      s2     <- sum(resids^2)/(N-Ki)                             # sigma squared
      if(is.null(R.restr)) {
        bcov   <- s2 * solve( crossprod( Xf ), tol=solvetol )     # coefficient covariance matrix
      } else {
        bcov   <- s2 * solve( W, tol=solvetol )[1:ncol(X),1:ncol(X)]  # coeff. covariance matrix
      }
    }
  }

  ## only for W2SLS estimation
  if(method=="W2SLS") {
    bl     <- b                       # coefficients of previous step
    bdif   <- b                       # difference of coefficients between this and previous step
    rcov   <- matrix( 0, G, G )                               # residual covariance matrix
    iter  <- 0
    while((sum(bdif^2)/sum(bl^2))^0.5>tol & iter < maxiter) {
      iter  <- iter+1
      bl    <- b                           # coefficients of previous step
      resids <- Y-X%*%b                     # residuals
      for(i in 1:G) residi[[i]] <- resids[(1+sum(n[1:i])-n[i]):(sum(n[1:i]))]
      for(i in 1:G)  {
        if(rcovformula==0) {
          rcov[i,i] <- sum(residi[[i]]*residi[[i]])/n[i]
        } else {
          rcov[i,i] <- sum(residi[[i]]*residi[[i]])/df[i]
        }
      }
      Oinv <- solve( rcov, tol=solvetol ) %x% diag(1,n[1],n[1]) # Omega inverse(= weight. matrix)
      if(is.null(R.restr)) {
        b <- solve(t(Xf) %*% Oinv %*% Xf, tol=solvetol) %*% t(Xf) %*% Oinv %*% Y  # (unrestr.) coeffic.
      } else {
        W <- rbind( cbind( t(Xf) %*% Oinv %*% Xf, t(R.restr) ),
                    cbind( R.restr, matrix(0, nrow(R.restr), nrow(R.restr))))
        V <- rbind( t(Xf) %*% Oinv %*% Y , q.restr )
        Winv <- solve( W, tol=solvetol )
        b <- ( Winv %*% V )[1:ncol(X)]     # restricted coefficients
      }
      bdif <- b - bl                       # difference of coefficients between this and previous step
    }
    if(is.null(R.restr)) {
      bcov <- solve(t(Xf) %*% Oinv %*% Xf, tol=solvetol )   # coefficient covariance matrix
    } else {
      bcov   <- Winv[1:ncol(X),1:ncol(X)]     # coefficient covariance matrix
    }
    resids <- Y - X %*% b                        # residuals
    for(i in 1:G) residi[[i]] <- resids[(1+sum(n[1:i])-n[i]):(sum(n[1:i]))]
  }

  ## only for 3SLS estimation
  if(method=="3SLS") {
    bl     <- b                       # coefficients of previous step
    bdif   <- b                       # difference of coefficients between this and previous step
    rcov   <- matrix( 0, G, G )                               # residual covariance matrix
    iter  <- 0
    while((sum(bdif^2)/sum(bl^2))^0.5>tol & iter < maxiter) {
      iter  <- iter+1
      bl    <- b                           # coefficients of previous step
      resids <- Y-X%*%b                     # residuals
      for(i in 1:G) residi[[i]] <- resids[(1+sum(n[1:i])-n[i]):(sum(n[1:i]))]
      for(i in 1:G) {
        for(j in 1:G) {
          if(rcovformula==0 | rcovformula==1) {
            rcov[i,j] <- sum(residi[[i]]*residi[[j]])/(
                          sqrt((n[i]-rcovformula*ki[i])*(n[j]-rcovformula*ki[j])))
          } else {
            rcov[i,j] <- sum(residi[[i]]*residi[[j]])/(
                          n[i]-ki[i]-ki[j] + sum( diag(
                          solve( crossprod( x[[i]] ), tol=solvetol ) %*%
                          crossprod( x[[i]], x[[j]]) %*%
                          solve( crossprod( x[[j]] ), tol=solvetol ) %*%
                          crossprod( x[[j]], x[[i]]))))
          }
        }
      }
      Oinv <- solve( rcov, tol=solvetol ) %x% diag(1,n[1],n[1])  # Omega inverse (= weighting matrix)
      if(formula3sls=="GLS") {
        if(is.null(R.restr)) {
          b <- solve(t(Xf) %*% Oinv %*% Xf, tol=solvetol) %*% t(Xf) %*% Oinv %*% Y  # (unrestr.) coeffic.
        } else {
          W <- rbind( cbind( t(Xf) %*% Oinv %*% Xf, t(R.restr) ),
                      cbind( R.restr, matrix(0, nrow(R.restr), nrow(R.restr))))
          V <- rbind( t(Xf) %*% Oinv %*% Y , q.restr )
          Winv <- solve( W, tol=solvetol )
          b <- ( Winv %*% V )[1:ncol(X)]     # restricted coefficients
        }
      }
      if(formula3sls=="IV") {
        if(is.null(R.restr)) {
          b <- solve(t(Xf) %*% Oinv %*% X, tol=solvetol) %*% t(Xf) %*% Oinv %*% Y  # (unrestr.) coeffic.
        } else {
          W <- rbind( cbind( t(Xf) %*% Oinv %*% X, t(R.restr) ),
                      cbind( R.restr, matrix(0, nrow(R.restr), nrow(R.restr))))
          V <- rbind( t(Xf) %*% Oinv %*% Y , q.restr )
          Winv <- solve( W, tol=solvetol )
          b <- ( Winv %*% V )[1:ncol(X)]     # restricted coefficients
        }
      }
      if(formula3sls=="GMM") {
        if(is.null(R.restr)) {
          b <- solve(t(X) %*% H %*% solve( t(H) %*% ( rcov %x% diag(1,n[1],n[1])) %*% H, tol=solvetol)
                  %*% t(H) %*% X, tol=solvetol) %*% t(X) %*% H %*% solve( t(H) %*%
                  ( rcov %x% diag(1,n[1],n[1])) %*% H, tol=solvetol) %*% t(H) %*% Y  #(unrestr.) coeffic.
        } else {
          W <- rbind( cbind( t(X) %*% H %*% solve( t(H) %*% ( rcov %x% diag(1,n[1],n[1]))
                              %*% H, tol=solvetol) %*% t(H) %*% X, t(R.restr) ),
                      cbind( R.restr, matrix(0, nrow(R.restr), nrow(R.restr))))
          V <- rbind( t(X) %*% H %*% solve( t(H) %*% ( rcov %x% diag(1,n[1],n[1]))
                      %*% H, tol=solvetol) %*% t(H) %*% Y , q.restr )
          Winv <- solve( W, tol=solvetol )
          b <- ( Winv %*% V )[1:ncol(X)]     # restricted coefficients
        }
      }
      if(formula3sls=="Schmidt") {
        if(is.null(R.restr)) {
          b <- solve( t(Xf) %*% Oinv %*% Xf, tol=solvetol) %*% ( t(Xf) %*% Oinv
                      %*% H %*% solve( crossprod( H ), tol=solvetol ) %*% crossprod(H, Y) )
                           # (unrestr.) coeffic.
        } else {
          W <- rbind( cbind( t(Xf) %*% Oinv %*% XF, t(R.restr) ),
                      cbind( R.restr, matrix(0, nrow(R.restr), nrow(R.restr))))
          V <- rbind( t(Xf) %*% Oinv %*% H %*% solve( crossprod( H ), tol=solvetol ) %*%
                      crossprod( H, Y ), q.restr )
          Winv <- solve( W, tol=solvetol )
          b <- ( Winv %*% V )[1:ncol(X)]     # restricted coefficients
        }
      }
      if(formula3sls=="EViews") {
        if(is.null(R.restr)) {
          b  <- b2 + solve(t(Xf) %*% Oinv %*% Xf, tol=solvetol) %*% ( t(Xf) %*% Oinv
                    %*% (Y -  X %*% b2) )   # (unrestr.) coeffic.
        } else {
          W <- rbind( cbind( t(Xf) %*% Oinv %*% Xf, t(R.restr) ),
                      cbind( R.restr, matrix(0, nrow(R.restr), nrow(R.restr))))
          V <- rbind( t(Xf) %*% Oinv %*% (Y -  X %*% b2) , q.restr )
          Winv <- solve( W, tol=solvetol )
          b <- b2 + ( Winv %*% V )[1:ncol(X)]     # restricted coefficients
        }
      }
      bdif <- b - bl                       # difference of coefficients between this and previous step
    }
    if(formula3sls=="GLS") {
      if(is.null(R.restr)) {
        bcov <- solve(t(Xf) %*% Oinv %*% Xf, tol=solvetol )   # coefficient covariance matrix
      } else {
        bcov   <- Winv[1:ncol(X),1:ncol(X)]     # coefficient covariance matrix
      }
    }
    if(formula3sls=="IV") {
      if(is.null(R.restr)) {
        bcov <- solve(t(Xf) %*% Oinv %*% X, tol=solvetol )   # final step coefficient covariance matrix
      } else {
        bcov   <- Winv[1:ncol(X),1:ncol(X)]     # coefficient covariance matrix
      }
    }
    if(formula3sls=="GMM") {
      if(is.null(R.restr)) {
        bcov <- solve(t(Xf) %*% Oinv %*% Xf, tol=solvetol )   # final step coefficient covariance matrix
      } else {
        bcov   <- Winv[1:ncol(X),1:ncol(X)]     # coefficient covariance matrix
      }
    }
    if(formula3sls=="Schmidt") {
      if(is.null(R.restr)) {
        bcov <- solve(t(Xf) %*% Oinv %*% Xf, tol=solvetol )   # final step coefficient covariance matrix
      } else {
        bcov   <- Winv[1:ncol(X),1:ncol(X)]     # coefficient covariance matrix
      }
    }
    if(formula3sls=="EViews") {
      if(is.null(R.restr)) {
        bcov <- solve(t(Xf) %*% Oinv %*% Xf, tol=solvetol )   # final step coefficient covariance matrix
      } else {
        W <- rbind( cbind( t(Xf) %*% Oinv %*% Xf, t(R.restr) ), cbind( R.restr, matrix(0,K-Ki,K-Ki)))
        V <- rbind( t(Xf) %*% Oinv %*% Y , q.restr )
        bcov <- solve( W, tol=solvetol )[1:ncol(X),1:ncol(X)]     # coefficient covariance matrix
      }
    }
    resids <- Y - X %*% b                        # residuals
    for(i in 1:G) residi[[i]] <- resids[(1+sum(n[1:i])-n[i]):(sum(n[1:i]))]
  }

  ## for all estimation methods
  pred  <- X %*% b                              # predicted endogenous values
  bt    <- NULL
  btcov <- NULL
  if(!is.null(TX)) {
    bt <- b
    b  <- TX %*% bt
    btcov <- bcov
    bcov  <- TX %*% btcov %*% t(TX)
  }
  se     <- diag(bcov)^0.5                       # standard errors of all estimated coefficients
  t      <- b/se                                 # t-values of all estimated coefficients
  if(probdfsys) {
    prob <- 2*( 1-pt(abs(t), N - Ki))            # p-values of all estimated coefficients
  } else {
    prob <- matrix( 0, 0, 1 )                    # p-values of all estimated coefficients
  }



  ## equation wise results
  for(i in 1:G) {
    bi     <- b[(1+sum(k[1:i])-k[i]):(sum(k[1:i]))]       # estimated coefficients of equation i
    sei    <- c(se[(1+sum(k[1:i])-k[i]):(sum(k[1:i]))])   # std. errors of est. param. of equation i
    ti     <- c(t[(1+sum(k[1:i])-k[i]):(sum(k[1:i]))])    # t-values of estim. param. of equation i
    bcovi  <- bcov[(1+sum(k[1:i])-k[i]):(sum(k[1:i])),(1+sum(k[1:i])-k[i]):(sum(k[1:i]))]
                                        # covariance matrix of estimated coefficients of equation i
    bi     <- array(bi,c(k[i],1))
    rownames(bi) <- colnames(x[[i]])
    attr(bi,"names") <- colnames(x[[i]])

    if(probdfsys) {
      probi <- c(prob[(1+sum(k[1:i])-k[i]):(sum(k[1:i]))]) # p-values of estim. param. of equation i
    } else {
      probi <- 2*( 1 - pt(abs(ti), df[i] ))              # p-values of estim. param. of equation i
      prob <- c(prob,probi)                              # p-values of all estimated coefficients
    }
    ssr    <- sum(residi[[i]]^2)                         # sum of squared residuals
    mse    <- ssr/df[i]                                  # estimated variance of residuals
    rmse   <- sqrt( mse )                                # estimated standard error of residuals
    r2     <- 1 - ssr/(t(y[[i]])%*%y[[i]]-n[i]*mean(y[[i]])^2)
    adjr2  <- 1 - ((n[i]-1)/df[i])*(1-r2)
    predi  <- pred[(1+sum(n[1:i])-n[i]):(sum(n[1:i]))]
    datai  <- model.frame( eqns[[i]] )
    if(method=="2SLS" | method=="3SLS") {
      datai <- cbind( datai, model.frame( instl[[i]] ))
    }

    if(i==1) {
      alldata <- datai                    # dataframe for all data used for estimation
    } else {
      alldata <- cbind( alldata, datai )  # dataframe for all data used for estimation
    }

    ## build the "return" structure for the equations
    resulti$method       <- method
    resulti$i            <- i               # equation number
    resulti$eqnlabel     <- eqnlabels[[i]]
    resulti$formula      <- eqns[[i]]
    resulti$n            <- n[i]            # number of observations
    resulti$k            <- k[i]            # number of coefficients/regressors
    resulti$ki           <- ki[i]           # number of linear independent coefficients
    resulti$df           <- df[i]           # degrees of freedom of residuals
    resulti$b            <- c(bi)           # estimated coefficients
    resulti$se           <- c(sei)          # standard errors of estimated coefficients
    resulti$t            <- c(ti)           # t-values of estimated coefficients
    resulti$p            <- c(probi)        # p-values of estimated coefficients
    resulti$bcov         <- bcovi           # covariance matrix of estimated coefficients
    resulti$y            <- y[[i]]          # vector of endogenous variables
    resulti$x            <- x[[i]]          # matrix of regressors
    resulti$data         <- datai           # data frame of this equation (incl. instruments)
    resulti$predicted    <- predi           # predicted values
    resulti$residuals    <- residi[[i]]     # residuals
    resulti$ssr          <- ssr             # sum of squared errors/residuals
    resulti$mse          <- mse             # estimated variance of the residuals (mean squared error)
    resulti$s2           <- mse             #        the same (sigma hat squared)
    resulti$rmse         <- rmse            # estimated standard error of the residuals
    resulti$s            <- rmse            #        the same (sigma hat)
    resulti$r2           <- r2              # R-sqared value
    resulti$adjr2        <- adjr2           # adjusted R-squared value
    if(method=="2SLS" | method=="3SLS") {
      resulti$inst         <- instl[[i]]
      resulti$h            <- h[[i]]          # matrix of instrumental variables
    }
    class(resulti)	     <- "systemfit.equation"
    results$eq[[i]]      <- resulti
  }

  ## results of the total system
  olsr2 <- 1 - t(resids) %*% resids / ( t(Y) %*% ( diag(1,G,G)     # OLS system R2
               %x% ( diag( 1, n[1], n[1]) - rep(1, n[1]) %*% t(rep(1,n[1])) / n[1])) %*% Y)
  if(method=="SUR" | method=="3SLS") {
    rcovest <- rcov                   # residual covariance matrix used for estimation
  }
  rcov    <- matrix( 0, G, G )                        # residual covariance matrix
  rcor    <- matrix( 0, G, G )                        # residual covariance matrix
  for(i in 1:G) { for(j in 1:G) {
      rcor[i,j] <- sum( residi[[i]] * residi[[j]] ) / (
                     sqrt( sum( residi[[i]]^2 ) * sum( residi[[j]]^2 )))
      if(rcovformula==0 | rcovformula==1) {
        rcov[i,j] <- sum(residi[[i]]*residi[[j]])/(
                     sqrt((n[i]-rcovformula*ki[i])*(n[j]-rcovformula*ki[j])))
      } else {
        rcov[i,j] <- sum(residi[[i]]*residi[[j]])/(
                           n[i]-ki[i]-ki[j] + sum( diag(
                           solve( crossprod( x[[i]] ), tol=solvetol ) %*%
                           crossprod( x[[i]], x[[j]] ) %*%
                           solve( crossprod( x[[j]] ), tol=solvetol ) %*%
                           crossprod( x[[j]], x[[i]] ))))
             }

  } }
  drcov <- det(rcov, tol=solvetol)
   mcelr2 <- 1 - ( t(resids) %*% ( solve(rcov, tol=solvetol) %x% diag(1, n[1],n[1])) %*% resids ) / (
              t(Y) %*% ( solve(rcov, tol=solvetol) %x% ( diag(1,n[1],n[1]) - rep(1,n[1]) %*%
              t(rep(1,n[1])) / n[1] )) %*% Y )   # McElroy's (1977a) R2

  b              <- array(b,c(K,1))
  rownames(b)    <- xnames
  colnames(b)    <- c("coefficient")
  se             <- array(se,c(K,1))
  rownames(se)   <- xnames
  colnames(se)   <- c("standard error")
  t              <- array(t,c(K,1))
  rownames(t)    <- xnames
  colnames(t)    <- c("t-statistic")
  prob           <- array(prob,c(K,1))
  rownames(prob) <- xnames
  colnames(prob) <- c("p-value")


  ## build the "return" structure for the whole system
  results$method  <- method
  results$g       <- G              # number of equations
  results$n       <- N              # total number of observations
  results$k       <- K              # total number of coefficients
  results$ki      <- Ki             # total number of linear independent coefficients
  results$b       <- b              # all estimated coefficients
  results$bt      <- bt             # transformed vector of estimated coefficients
  results$se      <- se             # standard errors of estimated coefficients
  results$t       <- t              # t-values of estimated coefficients
  results$p       <- prob           # p-values of estimated coefficients
  results$bcov    <- bcov           # coefficients covariance matrix
  results$btcov   <- btcov          # covariance matrix for transformed coefficient vector
  results$rcov    <- rcov           # residual covarance matrix
  results$drcov   <- drcov          # determinant of residual covarance matrix
  results$rcor    <- rcor           # residual correlation matrix
  results$olsr2   <- olsr2          # R-squared value of the equation system
  results$iter    <- iter           # residual correlation matrix
  results$y       <- y              # vector of all (stacked) endogenous variables
  results$x       <- X              # matrix of all (diagonally stacked) regressors
  results$resids  <- resids         # vector of all (stacked) residuals
  results$data    <- alldata        # data frame for all data used in the system estimation
  if(method=="2SLS" | method=="3SLS") {
    results$h       <- H            # matrix of all (diagonally stacked) instrumental variables
  }
  if(method=="SUR" | method=="3SLS") {
    results$rcovest <- rcovest      # residual covarance matrix used for estimation
    results$mcelr2  <- mcelr2       # McElroy's R-squared value for the equation system
  }
  results$R.restr <- R.restr
  results$q.restr <- q.restr
  results$TX      <- TX
  results$maxiter <- maxiter
  results$tol     <- tol
  results$rcovformula     <- rcovformula
  results$formula3sls     <- formula3sls
  results$probdfsys       <- probdfsys
  results$single.eq.sigma <- single.eq.sigma
  results$solvetol        <- solvetol
  class(results)  <- "systemfit.system"

  detach(data)
  results
}


## print the (summary) results that belong to the whole system
summary.systemfit.system <- function(object,...) {
  summary.systemfit.system <- object
  summary.systemfit.system
}

## print the results that belong to the whole system
print.systemfit.system <- function( x, digits=6,... ) {
  object <- x

  save.digits <- unlist(options(digits=digits))
  on.exit(options(digits=save.digits))

  table <- NULL
  labels <- NULL

  cat("\n")
  cat("systemfit results \n")
  cat("method: ")
  if(!is.null(object$iter)) if(object$iter>1) cat("iterated ")
  cat( paste( object$method, "\n\n"))
  if(!is.null(object$iter)) {
    if(object$iter>1) {
      if(object$iter<object$maxiter) {
        cat( paste( "convergence achieved after",object$iter,"iterations\n\n" ) )
      } else {
        cat( paste( "warning: convergence not achieved after",object$iter,"iterations\n\n" ) )
      }
    }
  }
  for(i in 1:object$g) {
    row <- NULL
    row <- cbind( round( object$eq[[i]]$n,     digits ),
                  round( object$eq[[i]]$df,    digits ),
                  round( object$eq[[i]]$ssr,   digits ),
                  round( object$eq[[i]]$mse,   digits ),
                  round( object$eq[[i]]$rmse,  digits ),
                  round( object$eq[[i]]$r2,    digits ),
                  round( object$eq[[i]]$adjr2, digits ))
    table  <- rbind( table, row )
    labels <- rbind( labels, object$eq[[i]]$eqnlabel )
  }
  rownames(table) <- c( labels )
  colnames(table) <- c("N","DF", "SSR", "MSE", "RMSE", "R2", "Adj R2" )

  print.matrix(table, quote = FALSE, right = TRUE )
  cat("\n")

  if(!is.null(object$rcovest)) {
    cat("The covariance matrix of the residuals used for estimation\n")
    rcov <- object$rcovest
    rownames(rcov) <- labels
    colnames(rcov) <- labels
    print( rcov )
    cat("\n")
    if( min(eigen( object$rcov, only.values=TRUE)$values) < 0 ) {
      cat("warning: this covariance matrix is NOT positive semidefinit!\n")
      cat("\n")
    }
  }

  cat("The covariance matrix of the residuals\n")
  rcov <- object$rcov
  rownames(rcov) <- labels
  colnames(rcov) <- labels
  print( rcov )
  cat("\n")

  cat("The correlations of the residuals\n")
  rcor <- object$rcor
  rownames(rcor) <- labels
  colnames(rcor) <- labels
  print( rcor )
  cat("\n")

  cat("The determinant of the residual covariance matrix: ")
  cat(object$drcov)
  cat("\n")

  cat("OLS R-squared value of the system: ")
  cat(object$olsr2)
  cat("\n")

  if(!is.null(object$mcelr2)) {
    cat("McElroy's R-squared value for the system: ")
    cat(object$mcelr2)
    cat("\n")
  }
  ## now print the individual equations
  for(i in 1:object$g) {
      print( object$eq[[i]], digits )
  }
}


## print the (summary) results for a single equation
summary.systemfit.equation <- function(object,...) {
  summary.systemfit.equation <- object
  summary.systemfit.equation
}


## print the results for a single equation
print.systemfit.equation <- function( x, digits=6, ... ) {
  object <- x

  save.digits <- unlist(options(digits=digits))
  on.exit(options(digits=save.digits))

  cat("\n")
  cat( paste( object$method, "estimates for", object$eqnlabel, " (equation", object$i, ")\n" ) )

  cat("Model Formula: ")
  print(object$formula)
  if(!is.null(object$inst)) {
    cat("Instruments: ")
    print(object$inst)
  }
  cat("\n")

  Signif <- symnum(object$p, corr = FALSE, na = FALSE,
                   cutpoints = c(0,  .001,.01,.05, .1, 1),
                   symbols   = c("***","**","*","."," "))

  table <- cbind(round( object$b,  digits ),
                 round( object$se, digits ),
                 round( object$t,  digits ),
                 round( object$p,  digits ),
                 Signif)

  rownames(table) <- names(object$b)
  colnames(table) <- c("Estimate","Std. Error","t value","Pr(>|t|)","")

  print.matrix(table, quote = FALSE, right = TRUE )
  cat("---\nSignif. codes: ",attr(Signif,"legend"),"\n")

  cat(paste("\nResidual standard error:", round(object$s, digits),  ## s ist the variance, isn't it???
            "on", object$df, "degrees of freedom\n"))

  cat( paste( "Number of observations:", round(object$n, digits),
              "Degrees of Freedom:", round(object$df, digits),"\n" ) )

  cat( paste( "SSR:",     round(object$ssr,    digits),
              "MSE:", round(object$mse, digits),
              "Root MSE:",   round(object$rmse,  digits), "\n" ) )

  cat( paste( "Multiple R-Squared:", round(object$r2,    digits),
              "Adjusted R-Squared:", round(object$adjr2, digits),
              "\n" ) )
  cat("\n")
}


## calculate predicted values and its standard errors and limits
prediction.systemfit <- function( object, data=object$data, alpha=0.05) {
   attach(data)
   results <- list()
   g       <- object$g
   n       <- array(NA,c(g))
   eqns    <- list()
   x       <- list()               # regressors equation-wise
   X       <- matrix( 0, 0, 0 )    # stacked matrices of all regressors (unrestricted)
   for(i in 1:g )  {
      eqns[[i]] <- object$eq[[i]]$formula
      x[[i]] <-  model.matrix( eqns[[i]] )
      X      <-  rbind( cbind( X, matrix( 0, nrow( X ), ncol( x[[i]] ))),
                       cbind( matrix( 0, nrow( x[[i]] ), ncol( X )), x[[i]]))
      n[i]   <-  nrow( x[[i]] )
   }
   Y <- X %*% object$b
   if(object$method=="SUR" | object$method=="3SLS") {
      ycov <- X %*% object$bcov %*% t(X) + object$rcov %x% diag(1,n[1],n[1])
   }
   for(i in 1:g) {
      resulti <- list()
      Yi <- Y[(1+sum(n[1:i])-n[i]):sum(n[1:i]),]
      if(object$method=="SUR" | object$method=="3SLS") {
         ycovi <- ycov[(1+sum(n[1:i])-n[i]):sum(n[1:i]),(1+sum(n[1:i])-n[i]):sum(n[1:i])]
      } else {
         ycovi <- x[[i]] %*% object$eq[[i]]$bcov %*% t(x[[i]]) + object$eq[[i]]$s2
      }
      if(nrow(data)==1) {
         sei    <- sqrt( ycovi )
      } else {
         sei    <- sqrt( diag(ycovi) )
      }
      tval   <- qt( 1 - alpha/2, object$eq[[i]]$df )
      limits <- cbind( Yi - (tval*sei), Yi + (tval*sei) )

      resulti$predicted          <- Yi
      resulti$se.prediction      <- sei
      resulti$prediction.limits  <- limits
      results[[i]]               <- resulti
   }
   detach(data)
   results
}


## this function returns a vector of the
## cross-equation corrlations between eq i and eq j
## from the results set for equation ij
correlation.systemfit <- function( results, eqni, eqnj ) {
  k <- array( 0, c(results$g))
  for(i in 1:results$g) k[i] <- results$eq[[i]]$k
  cij <- results$bcov[(1+sum(k[1:eqni])-k[eqni]):(sum(k[1:eqni])),
                      (1+sum(k[1:eqnj])-k[eqnj]):(sum(k[1:eqnj]))]
  cii <- results$bcov[(1+sum(k[1:eqni])-k[eqni]):(sum(k[1:eqni])),
                      (1+sum(k[1:eqni])-k[eqni]):(sum(k[1:eqni]))]
  cjj <- results$bcov[(1+sum(k[1:eqnj])-k[eqnj]):(sum(k[1:eqnj])),
                      (1+sum(k[1:eqnj])-k[eqnj]):(sum(k[1:eqnj]))]
  rij <- NULL

  for(i in 1:results$eq[[1]]$n ) {
    xik    <- results$eq[[eqni]]$x[i,]
    xjk    <- results$eq[[eqnj]]$x[i,]
    top    <- xik %*% cij %*% xjk
    bottom <- sqrt( ( xik %*% cii %*% xik ) * ( xjk %*% cjj %*% xjk ) )
    rijk   <- top / bottom
    rij    <- rbind( rij, rijk )
  }
  rij
}

## determines the improvement of resultsj (3sls) over
## resultsi (2sls) for equation i and returns a matrix
## of the values, so you can examine the range, mean, etc
se.ratio.systemfit <- function( resultsi, resultsj, eqni ) {
  ratio <- NULL
  for(i in 1:resultsi$eq[[eqni]]$n ) {
    xik    <- resultsi$eq[[eqni]]$x[i,]
    top    <- sqrt( xik %*% resultsi$eq[[eqni]]$bcov %*% xik )
    bottom <- sqrt( xik %*% resultsj$eq[[eqni]]$bcov %*% xik )
    rk     <- top / bottom
    ratio  <- rbind( ratio, rk )
  }
  ratio
}


## this function returns test statistic for
## the hausman test which.... i forget, but people want to see it...
## from the sas docs
## given 2 estimators, b0 abd b1, where under the null hypothesis,
## both are consistent, but only b0 is asympt. efficient and
## under the alter. hypo only b1 is consistent, so the statistic (m) is

## man is this wrong...
hausman.systemfit <- function( li.results, fi.results ) {

  ## build the variance-covariance matrix
  ## for the full information and the limited information
  ## matricies
  ficovb <- NULL
  licovb <- NULL
  lib    <- NULL
  fib    <- NULL

  ## build the final large matrix...
  for(i in 1:li.results$g ) {
    fitr <- NULL
    litr <- NULL

    ## get the dimensions of the current matrix
    for(j in 1:li.results$g ) {
      if( i == j ) {
        litr <- cbind( litr, li.results$eq[[i]]$bcov )
      } else {
        ## bind the zero matrix to the row
        di   <- dim( li.results$eq[[i]]$bcov )[1]
        dj   <- dim( li.results$eq[[j]]$bcov )[1]
        litr <- cbind( litr, matrix( 0, di, dj ) )
      }
    }
    licovb <- rbind( licovb, litr )

    ## now add the rows of the parameter estimates
    ## to the big_beta matrix to compute the differences
    lib <- c( lib, li.results$eq[[i]]$b )
    fib <- c( fib, fi.results$eq[[i]]$b )
  }
  vli <- licovb
  vfi <- fi.results$bcov
  q   <- fib - lib

  hausman <- t( q ) %*% solve( vli - vfi ) %*% q
  hausman
}


## Likelihood Ratio Test
lrtest.systemfit <- function( resultc, resultu ) {
  lrtest <- list()
  if(resultc$method=="SUR" & resultu$method=="SUR") {
    n   <- resultu$eq[[1]]$n
    lrtest$df  <- resultu$ki - resultc$ki
    if(resultc$rcovformula != resultu$rcovformula) {
      stop("both estimations must use the same formula to calculate the residual covariance matrix!")
    }
    if(resultc$rcovformula == 0) {
      lrtest$lr  <- n * ( log( resultc$drcov ) - log( resultu$drcov ) )
    } else {
      residc <- array(resultc$resids,c(n,resultc$g))
      residu <- array(resultu$resids,c(n,resultu$g))
      lrtest$lr <- n * ( log( det( (t(residc) %*% residc)) ) - log( det( (t(residu) %*% residu))))
    }
    lrtest$p <- 1-pchisq( lrtest$lr, lrtest$df )
  }
  lrtest
}


