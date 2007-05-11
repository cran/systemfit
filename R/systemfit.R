###   $Id: systemfit.R 354 2007-03-28 14:48:14Z henningsena $
###
###            Simultaneous Equation Estimation for R
###
### Copyright 2002-2004 Jeff D. Hamann <jeff.hamann@forestinformatics.com>
###                     Arne Henningsen <http://www.arne-henningsen.de>
###
### This file is part of the nlsystemfit library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA


systemfit <- function( method,
                        eqns,
                        eqnlabels=names(eqns),
                        inst=NULL,
                        data=list(),
                        R.restr=NULL,
                        q.restr=matrix(0,max(nrow(R.restr),0),1),
                        TX=NULL,
                        maxiter=1,
                        tol=1e-5,
                        rcovformula=1,
                        centerResiduals = FALSE,
                        formula3sls="GLS",
                        probdfsys=!(is.null(R.restr) & is.null(TX)),
                        single.eq.sigma=(is.null(R.restr) & is.null(TX)),
                        solvetol=.Machine$double.eps,
                        saveMemory = ( nrow( data ) * length( eqns ) > 1000 &&
                           length( data ) > 0 ) )
{

   ## some tests
   if(!( method=="OLS" | method=="WLS" | method=="SUR" | method=="WSUR" |
         method=="2SLS" | method=="W2SLS" | method=="3SLS" | method=="W3SLS" )){
      stop( "The method must be 'OLS', 'WLS', 'SUR', 'WSUR',",
         " '2SLS', 'W2SLS', '3SLS', or 'W3SLS'" )
   }
   if( ( method=="2SLS" | method=="W2SLS" | method=="3SLS" | method=="W3SLS" ) &
         is.null(inst) ) {
      stop( "The methods '2SLS', 'W2SLS', '3SLS', and 'W3SLS' need instruments!" )
   }

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
  k       <- array( 0, c(G) )     # number of (unrestricted) coefficients/
                                  # regressors in each equation
  instl   <- list()               # list of the instruments for each equation
  ssr     <- array( 0, c(G))      # sum of squared residuals of each equation
  mse     <- array( 0, c(G))      # mean square error (residuals) of each equation
  rmse    <- array( 0, c(G))      # root of mse
  r2      <- array( 0, c(G))      # R-squared value
  adjr2   <- array( 0, c(G))      # adjusted R-squared value
  xnames  <- NULL                 # names of regressors

#   for(i in 1:G )  {
#     y[[i]] <-  eval( attr( terms( eqns[[i]] ), "variables" )[[2]] )
#     Y      <-  c( Y, y[[i]] )
#     x[[i]] <-  model.matrix( eqns[[i]] )
#     X      <-  rbind( cbind( X, matrix( 0, nrow( X ), ncol( x[[i]] ))),
#                        cbind( matrix( 0, nrow( x[[i]] ), ncol( X )), x[[i]]))
#     n[i]   <-  length( y[[i]] )
#     k[i]   <-  ncol(x[[i]])
#     for(j in 1:k[i]) {
#       xnames <- c( xnames, paste("eq",as.character(i),colnames( x[[i]] )[j] ))
#     }
#   }

   # the previous lines are subtituted by the following,
   # because Ott Toomet reported that they might lead to
   # problems with special data sets. He suggested the
   # following lines, which are copied from the survreg
   # package
   # how were we called?
   call <- match.call() # get the original call
   m0 <- match.call( expand.dots = FALSE ) #-"- without ...-expansion
   pos <- which( names( m0 ) == "data" )
   if( length( pos ) == 1 ) {
      results$data.name <- as.character( m0[ pos ] )
   } else {
      results$data.name <- "unknown"
   }
   rm( pos )
   temp <- c("", "data", "weights", "subset", "na.action")
                  # arguments for model matrices
   m0 <- m0[match(temp, names(m0), nomatch = 0)]
            # positions of temp-arguments
   m0[[1]] <- as.name("model.frame")
            # find matrices for individual models
   for(i in 1:G ) {
      m <- m0
      Terms <- terms(eqns[[i]], data = data)
      m$formula <- Terms
      m <- eval(m, parent.frame())
      weights <- model.extract(m, "weights")
      y[[i]] <- model.extract(m, "response")
      x[[i]] <- model.matrix(Terms, m)
      Y <- c(Y,y[[i]])
      X <- rbind( cbind( X, matrix( 0, nrow( X ), ncol( x[[i]] ))),
                  cbind( matrix( 0, nrow( x[[i]] ), ncol( X )), x[[i]]))
      n[i] <- length( y[[i]] )
      k[i] <- ncol(x[[i]])
      for(j in 1:k[i]) {
         xnames <- c( xnames, paste("eq",as.character(i),colnames( x[[i]] )[j] ))
      }
   }
   if( G > 1 ) {
      if( var ( n ) != 0 ) {
         stop( "Systems with unequal numbers of observations are not supported yet." )
      }
   }
   if( sum( n ) > 1000 && !saveMemory ) {
      warning( paste( "You have more than 1000 observations.",
         "Setting argument 'saveMemory' to TRUE speeds up",
         "the estimation. Estimation of larger data sets might even",
         "require this setting.\n" ) )
   }

   N  <- sum( n )    # total number of observations
   K  <- sum( k )    # total number of (unrestricted) coefficients/regressors
   Ki <- K           # total number of linear independent coefficients
   ki <- k           # total number of linear independent coefficients in each equation
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
            for(i in 1:G) {  # search for restrictions that are NOT cross-equation
               if( sum( R.restr[ j, (1+sum(k[1:i])-k[i]):(sum(k[1:i]))]^2) ==
                   sum(R.restr[j,]^2)) {
                  ki[i] <- ki[i]-1
               }
            }
         }
      }
    }
    df <- n - ki    # degress of freedom of each equation

  ## only for OLS, WLS and SUR estimation
  if(method=="OLS" | method=="WLS" | method=="SUR" | method=="WSUR") {
    if(is.null(R.restr)) {
      b <- solve( crossprod( X ), crossprod( X, Y ), tol=solvetol )
               # estimated coefficients
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
    if(single.eq.sigma) {
      rcov <- .calcRCov( resids, rcovformula = rcovformula, nObsEq = n,
         nCoefEq = ki, xEq = x, diag = TRUE, centered = centerResiduals,
         solvetol = solvetol )               # residual covariance matrix
      bcov <- .calcGLS( x = X, R.restr = R.restr, q.restr = q.restr, 
         sigma = rcov, nObsEq = n, solvetol = solvetol )
                    # coefficient covariance matrix
    } else {
      s2 <- .calcSigma2( resids, nObs = N, nCoef = Ki, rcovformula = rcovformula )
                           # sigma squared
      if(is.null(R.restr)) {
        bcov   <- s2 * solve( crossprod( X ), tol=solvetol )
                          # coefficient covariance matrix
      } else {
        bcov   <- s2 * solve( W, tol=solvetol )[1:ncol(X),1:ncol(X)]
                    # coefficient covariance matrix
      }
    }
  }

  ## only for WLS estimation
  if( method == "WLS" | method == "WSUR" ) {
    bl    <- b   # coefficients of previous step
    bdif  <- b   # difference of coefficients between this and previous step
    iter  <- 0
    while((sum(bdif^2)/sum(bl^2))^0.5>tol & iter < maxiter^( method == "WLS" ) ) {
      iter  <- iter+1
      bl    <- b                # coefficients of previous step
      resids <- Y - X %*% b     # residuals
      rcov <- .calcRCov( resids, rcovformula = rcovformula, nObsEq = n,
         nCoefEq = ki, xEq = x, diag = TRUE, centered = centerResiduals,
         solvetol = solvetol )
      b  <- .calcGLS( x = X, y = Y, R.restr = R.restr, q.restr = q.restr,
         sigma = rcov, nObsEq = n, solvetol = solvetol ) # coefficients
      bdif <- b-bl # difference of coefficients between this and previous step
    }
    bcov <- .calcGLS( x = X, R.restr = R.restr, q.restr = q.restr, 
       sigma = rcov, nObsEq = n, solvetol = solvetol )
       # final step coefficient covariance matrix
    resids <- Y - X %*% b                        # residuals
    for(i in 1:G) residi[[i]] <- resids[(1+sum(n[1:i])-n[i]):(sum(n[1:i]))]
  }

  ## only for SUR estimation
  if( method == "SUR" | method == "WSUR" ) {
    bl    <- b    # coefficients of previous step
    bdif  <- b    # difference of coefficients between this and previous step
    iter  <- 0
    while((sum(bdif^2)/sum(bl^2))^0.5>tol & iter < maxiter) {
      iter  <- iter+1
      bl    <- b                           # coefficients of previous step
      resids <- Y-X%*%b                     # residuals
      rcov <- .calcRCov( resids, rcovformula = rcovformula, nObsEq = n,
         nCoefEq = ki, xEq = x, centered = centerResiduals,
         solvetol = solvetol )
      b <- .calcGLS( x = X, y = Y, R.restr = R.restr, q.restr = q.restr, 
         sigma = rcov, nObsEq = n, solvetol = solvetol )     # coefficients
      bdif <- b-bl  # difference of coefficients between this and previous step
    }
    bcov <- .calcGLS( x = X, R.restr = R.restr, q.restr = q.restr, 
       sigma = rcov, nObsEq = n, solvetol = solvetol )
            # final step coefficient covariance matrix
    resids <- Y - X %*% b                        # residuals
    for(i in 1:G) residi[[i]] <- resids[(1+sum(n[1:i])-n[i]):(sum(n[1:i]))]
  }

  ## only for 2SLS, W2SLS and 3SLS estimation
  if( method == "2SLS" | method == "W2SLS" | method == "3SLS" |
      method == "W3SLS" ) {
    for(i in 1:G) {
      if(is.list(inst)) {
         instl[[i]] <- inst[[i]]
      } else {
         instl[[i]] <- inst
      }
    }
    Xf <- array(0,c(0,ncol(X)))       # fitted X values
    H  <- matrix( 0, 0, 0 )           # stacked matrices of all instruments
    h  <- list()
    for(i in 1:G) {
      Xi <- X[(1+sum(n[1:i])-n[i]):(sum(n[1:i])),]
            # regressors of the ith equation (including zeros)
      #h[[i]] <- model.matrix( instl[[i]] )
      # the following lines have been substituted for the previous
      # line due to changes in the data handling.
      # code provided by Ott Toomet
      m <- m0
      Terms <- terms(instl[[i]], data = data)
      m$formula <- Terms
      m <- eval(m, parent.frame())
      h[[i]] <- model.matrix(Terms, m)
      if( nrow( h[[ i ]] ) != nrow( Xi ) ) {
         stop( paste( "The instruments and the regressors of equation",
            as.character( i ), "have different numbers of observations." ) )
      }
      # extract instrument matrix
      Xf <- rbind(Xf, h[[i]] %*% solve( crossprod( h[[i]]) , tol=solvetol )
              %*% crossprod( h[[i]], Xi ))       # 'fitted' X-values
      H  <-  rbind( cbind( H, matrix( 0, nrow( H ), ncol( h[[i]] ))),
                         cbind( matrix( 0, nrow( h[[i]] ), ncol( H )), h[[i]]))

    }
    if(is.null(R.restr)) {
      b <- solve( crossprod( Xf ), crossprod( Xf, Y ), tol=solvetol )
         # 2nd stage coefficients
    } else {
      W <- rbind( cbind( crossprod(Xf), t(R.restr) ),
                  cbind( R.restr, matrix(0, nrow(R.restr), nrow(R.restr))))
      V <- rbind( t(Xf) %*% Y , q.restr )
      b <- ( solve( W, tol=solvetol ) %*% V )[1:ncol(X)] # restricted coefficients
    }
    b2 <- b
  }

  ## only for 2SLS estimation
  if(method=="2SLS") {
    resids <- Y - X %*% b                        # residuals
    if(single.eq.sigma) {
      rcov <- .calcRCov( resids, rcovformula = rcovformula, nObsEq = n,
         nCoefEq = ki, xEq = x, diag = TRUE, centered = centerResiduals,
         solvetol = solvetol )
      bcov <- .calcGLS( x = Xf, R.restr = R.restr, q.restr = q.restr, 
         sigma = rcov, nObsEq = n, solvetol = solvetol )  # coefficient covariance matrix
    } else {
      s2 <- .calcSigma2( resids, nObs = N, nCoef = Ki, rcovformula = rcovformula )
                           # sigma squared
      if(is.null(R.restr)) {
        bcov   <- s2 * solve( crossprod( Xf ), tol=solvetol )
                  # coefficient covariance matrix
      } else {
        bcov   <- s2 * solve( W, tol=solvetol )[1:ncol(X),1:ncol(X)]
                    # coeff. covariance matrix
      }
    }
  }

  ## only for W2SLS estimation
  if( method == "W2SLS" | method == "W3SLS" ) {
    bl     <- b   # coefficients of previous step
    bdif   <- b   # difference of coefficients between this and previous step
    iter  <- 0
    while((sum(bdif^2)/sum(bl^2))^0.5>tol & iter < maxiter^( method == "W2LS" ) ) {
      iter  <- iter+1
      bl    <- b                           # coefficients of previous step
      resids <- Y-X%*%b                     # residuals
      rcov <- .calcRCov( resids, rcovformula = rcovformula, nObsEq = n,
         nCoefEq = ki, xEq = x, diag = TRUE, centered = centerResiduals,
         solvetol = solvetol )
      b <- .calcGLS( x = Xf, y = Y, R.restr = R.restr, q.restr = q.restr, 
         sigma = rcov, nObsEq = n, solvetol = solvetol )          # (unrestr.) coeffic.
      bdif <- b - bl # difference of coefficients between this and previous step
    }
    bcov <- .calcGLS( x = Xf, R.restr = R.restr, q.restr = q.restr, 
       sigma = rcov, nObsEq = n, solvetol = solvetol )  # coefficient covariance matrix
    resids <- Y - X %*% b                        # residuals
    for(i in 1:G) residi[[i]] <- resids[(1+sum(n[1:i])-n[i]):(sum(n[1:i]))]
  }

  ## only for 3SLS estimation
  if( method == "3SLS" | method == "W3SLS" ) {
    bl     <- b  # coefficients of previous step
    bdif   <- b  # difference of coefficients between this and previous step
    iter  <- 0
    while((sum(bdif^2)/sum(bl^2))^0.5>tol & iter < maxiter) {
      iter  <- iter+1
      bl    <- b                           # coefficients of previous step
      resids <- Y-X%*%b                     # residuals
      rcov <- .calcRCov( resids, rcovformula = rcovformula, nObsEq = n,
         nCoefEq = ki, xEq = x, centered = centerResiduals, solvetol = solvetol )
      if(formula3sls=="GLS") {
         b <- .calcGLS( x = Xf, y = Y, R.restr = R.restr, q.restr = q.restr, 
            sigma = rcov, nObsEq = n, solvetol = solvetol )  # (unrestr.) coeffic.
      }
      if(formula3sls=="IV") {
         b <- .calcGLS( x = Xf, x2 = X, y = Y, R.restr = R.restr, q.restr = q.restr, 
            sigma = rcov, nObsEq = n, solvetol = solvetol )   # (unrestr.) coeffic.
      }
      if(formula3sls=="GMM") {
        HtOmega <- .calcXtOmegaInv( x = H, sigma = rcov, nObsEq = n,
           invertSigma = FALSE )
        if(is.null(R.restr)) {
          b <- solve(t(X) %*% H %*% solve( HtOmega %*%
                 H, tol=solvetol) %*% t(H) %*% X, tol=solvetol) %*% t(X) %*% H %*%
                 solve( HtOmega %*%
                 H, tol=solvetol) %*% t(H) %*% Y  #(unrestr.) coeffic.
        } else {
          W <- rbind( cbind( t(X) %*% H %*% solve( HtOmega
                              %*% H, tol=solvetol) %*% t(H) %*% X, t(R.restr) ),
                      cbind( R.restr, matrix(0, nrow(R.restr), nrow(R.restr))))
          V <- rbind( t(X) %*% H %*% solve( HtOmega
                      %*% H, tol=solvetol) %*% t(H) %*% Y , q.restr )
          Winv <- solve( W, tol=solvetol )
          b <- ( Winv %*% V )[1:ncol(X)]     # restricted coefficients
        }
      }
      if(formula3sls=="Schmidt") {
        XftOmegaInv <- .calcXtOmegaInv( x = Xf, sigma = rcov, nObsEq = n,
           solvetol = solvetol )
        if(is.null(R.restr)) {
          b <- solve( t(Xf) %*% t( XftOmegaInv ), tol=solvetol) %*% ( XftOmegaInv
                      %*% H %*% solve( crossprod( H ), tol=solvetol ) %*% crossprod(H, Y) )
                           # (unrestr.) coeffic.
        } else {
          W <- rbind( cbind( t(Xf) %*% t( XftOmegaInv ), t(R.restr) ),
                      cbind( R.restr, matrix(0, nrow(R.restr), nrow(R.restr))))
          V <- rbind( XftOmegaInv %*% H %*% solve( crossprod( H ), tol=solvetol ) %*%
                      crossprod( H, Y ), q.restr )
          Winv <- solve( W, tol=solvetol )
          b <- ( Winv %*% V )[1:ncol(X)]     # restricted coefficients
        }
      }
      if(formula3sls=="EViews") {
         b <- b2 + .calcGLS( x = Xf, y = ( Y -  X %*% b2 ), 
            R.restr = R.restr, q.restr = q.restr, sigma = rcov,
            nObsEq = n, solvetol = solvetol )  # (unrestr.) coeffic.
      }
      bdif <- b - bl # difference of coefficients between this and previous step
    }
    if(formula3sls=="GLS") {
       bcov <- .calcGLS( x = Xf, R.restr = R.restr, q.restr = q.restr, 
          sigma = rcov, nObsEq = n, solvetol = solvetol )  # coefficient covariance matrix
    }
    if(formula3sls=="IV") {
       bcov <- .calcGLS( x = Xf, x2 = X, R.restr = R.restr, q.restr = q.restr, 
          sigma = rcov, nObsEq = n, solvetol = solvetol )   
    }
    if(formula3sls=="GMM") {
      if(is.null(R.restr)) {
        bcov <- solve( t(X) %*% H %*% solve( HtOmega %*% H, tol=solvetol ) %*%
           t(H) %*% X, tol=solvetol )
                # final step coefficient covariance matrix
      } else {
        bcov   <- Winv[1:ncol(X),1:ncol(X)] # coefficient covariance matrix
      }
    }
    if(formula3sls=="Schmidt") {
      XftOmegaInv <- .calcXtOmegaInv( x = Xf, sigma = rcov, nObsEq = n,
         solvetol = solvetol )
      PH <- H %*%  solve( t(H) %*% H, tol=solvetol ) %*% t(H)
      PHOmega <- .calcXtOmegaInv( x = t( PH ), sigma = rcov, nObsEq = n,
           invertSigma = FALSE )
      if(is.null(R.restr)) {
         bcov <- solve( XftOmegaInv %*% Xf, tol=solvetol ) %*%
            XftOmegaInv %*% PHOmega %*%
            PH %*% t( XftOmegaInv ) %*% solve( XftOmegaInv %*% Xf, tol=solvetol )
                  # final step coefficient covariance matrix
      } else {
         VV <- XftOmegaInv %*% PHOmega %*%
            PH %*% t( XftOmegaInv )
         VV <- rbind( cbind( VV, matrix( 0, nrow( VV ), nrow( R.restr ) ) ),
            matrix( 0, nrow( R.restr ), nrow( VV ) + nrow( R.restr ) ) )
         bcov <- ( Winv %*% VV %*% Winv )[ 1:ncol(X), 1:ncol(X) ]
                  # coefficient covariance matrix
      }
    }
    if(formula3sls=="EViews") {
       bcov <- .calcGLS( x = Xf, R.restr = R.restr, q.restr = q.restr, 
          sigma = rcov, nObsEq = n, solvetol = solvetol )  # final step coefficient covariance matrix
    }
    resids <- Y - X %*% b                        # residuals
  }

  ## for all estimation methods
  fitted <- X %*% b                              # fitted endogenous values
  bt     <- NULL
  btcov  <- NULL
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
    residi[[i]] <- resids[ ( 1 + sum(n[1:i]) -n[i] ):( sum(n[1:i]) ) ]
    bi     <- b[(1+sum(k[1:i])-k[i]):(sum(k[1:i]))]
              # estimated coefficients of equation i
    sei    <- c(se[(1+sum(k[1:i])-k[i]):(sum(k[1:i]))])
              # std. errors of est. param. of equation i
    ti     <- c(t[(1+sum(k[1:i])-k[i]):(sum(k[1:i]))])
              # t-values of estim. param. of equation i
    bcovi  <- bcov[(1+sum(k[1:i])-k[i]):(sum(k[1:i])),(1+sum(k[1:i])-k[i]):(sum(k[1:i]))]
              # covariance matrix of estimated coefficients of equation i
    if(probdfsys) {
      probi <- c(prob[(1+sum(k[1:i])-k[i]):(sum(k[1:i]))])
               # p-values of estim. param. of equation i
    } else {
      probi <- 2*( 1 - pt(abs(ti), df[i] ))
               # p-values of estim. param. of equation i
      prob <- c(prob,probi) # p-values of all estimated coefficients
    }
    names( bi ) <- names( sei ) <- colnames( x[[i]] )
    names( ti ) <- names( probi ) <- colnames( x[[i]] )
    rownames( bcovi ) <- colnames( bcovi ) <- colnames( x[[i]] )

    ssr    <- sum(residi[[i]]^2)                         # sum of squared residuals
    mse    <- ssr/df[i]                                  # estimated variance of residuals
    rmse   <- sqrt( mse )                                # estimated standard error of residuals
    r2     <- 1 - ssr/(t(y[[i]])%*%y[[i]]-n[i]*mean(y[[i]])^2)
    adjr2  <- 1 - ((n[i]-1)/df[i])*(1-r2)
    fittedi <- fitted[(1+sum(n[1:i])-n[i]):(sum(n[1:i]))]
    #datai  <- model.frame( eqns[[i]] )
    # the following lines have to be substituted for the previous
    # line due to changes in the data handling.
    # code provided by Ott Toomet
    m <- m0
    Terms <- terms( eqns[[i]], data = data)
    m$formula <- Terms
    m <- eval(m, parent.frame())
    # datai <- model.frame(Terms, m)
    # the following lines are substituted for the previous line to
    # allow transformed variables in the formulas (e.g. "log(x1)")
    # code provided by Mikko Pakkanen
    resp <- model.extract( m, "response" )
    # using model.matrix instead of model.frame, need to get the output
    # variable separately
    datai <- data.frame( cbind( resp, ( model.matrix( Terms, m ) )[ , -1 ] ) )
    # I guess there's a better way to extract the name of the output variable?
    names( datai )[1] <- as.character( terms( eqns[[ i ]] ) )[2]
    rm( resp )
    if( method == "2SLS" | method == "W2SLS" | method == "3SLS" |
        method == "W3SLS" ) {
      #datai <- cbind( datai, model.frame( instl[[i]] ))
      # the following lines have to be substituted for the previous
      # line due to changes in the data handling.
      # code provided by Ott Toomet
      m <- m0
      Terms <- terms(instl[[i]], data = data)
      m$formula <- Terms
      m <- eval(m, parent.frame())
      # datai <- cbind( datai, model.frame(Terms, m))
      # the following lines are substituted for the previous line to
      # allow transformed variables in the formulas (e.g. "log(x1)")
      # code provided by Mikko Pakkanen
      datai <- cbind( datai,
         as.data.frame( ( model.matrix( Terms, m ) )[ , -1 ] ) )
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
    resulti$dfSys        <- N- Ki           # degrees of freedom of residuals of the whole system
    resulti$probdfsys    <- probdfsys       #
    resulti$b            <- c( bi )         # estimated coefficients
    resulti$se           <- c( sei )        # standard errors of estimated coefficients
    resulti$t            <- c( ti )         # t-values of estimated coefficients
    resulti$p            <- c( probi )      # p-values of estimated coefficients
    resulti$bcov         <- bcovi           # covariance matrix of estimated coefficients
    resulti$y            <- y[[i]]          # vector of endogenous variables
    resulti$x            <- x[[i]]          # matrix of regressors
    resulti$data         <- datai           # data frame of this equation (incl. instruments)
    resulti$fitted       <- fittedi         # fitted values
    resulti$residuals    <- residi[[i]]     # residuals
    resulti$ssr          <- ssr             # sum of squared errors/residuals
    resulti$mse          <- mse             # estimated variance of the residuals (mean squared error)
    resulti$s2           <- mse             #        the same (sigma hat squared)
    resulti$rmse         <- rmse            # estimated standard error of the residuals
    resulti$s            <- rmse            #        the same (sigma hat)
    resulti$r2           <- r2              # R-sqared value
    resulti$adjr2        <- adjr2           # adjusted R-squared value
    if( method == "2SLS" | method == "W2SLS" | method == "3SLS" |
        method == "W3SLS" ) {
      resulti$inst         <- instl[[i]]
      resulti$h            <- h[[i]]          # matrix of instrumental variables
    }
    class(resulti)        <- "systemfit.equation"
    results$eq[[i]]      <- resulti
  }

  ## results of the total system
  #olsr2 <- 1 - t(resids) %*% resids / ( t(Y) %*% ( diag(1,G,G)     # OLS system R2
  #             %x% ( diag( 1, n[1], n[1]) - rep(1, n[1]) %*% t(rep(1,n[1])) / n[1])) %*% Y)
  # the following lines are substituted for the previous 2 lines to increase
  # speed ( idea suggested by Ott Toomet )
   meanY <- numeric(length(Y)) # compute mean of Y by equations
   for(i in 1:G) {
      meanY[ (1+sum(n[1:i])-n[i]):(sum(n[1:i])) ] <-
         mean( Y[ (1+sum(n[1:i])-n[i]):(sum(n[1:i])) ])
   }
   olsr2 <- 1 - t(resids) %*% resids / sum( ( Y - meanY )^2 )
                        # OLS system R2
  if( method == "SUR" | method == "WSUR" | method == "3SLS" |
      method == "W3SLS" ) {
    rcovest <- rcov                   # residual covariance matrix used for estimation
  }
  rcov <- .calcRCov( resids, rcovformula = rcovformula, nObsEq = n,
      nCoefEq = ki, xEq = x, centered = centerResiduals, solvetol = solvetol )
  drcov <- det(rcov, tol=solvetol)
  if( !saveMemory ) {
#       # original formula from McElroy (1977)
#       mcelr2 <- 1 - ( t(resids) %*% ( solve(rcov, tol=solvetol) %x%
#                 diag(1, n[1],n[1])) %*% resids ) /
#                 ( t(Y) %*% ( solve(rcov, tol=solvetol ) %x%
#                 ( diag(1,n[1],n[1] ) - rep(1,n[1]) %*%
#                 t(rep(1,n[1])) / n[1] )) %*% Y )   # McElroy's (1977a) R2
      # first formula from Greene (2003, p. 345) (numerator modified to save memory)
      rtOmega <- .calcXtOmegaInv( x = resids, sigma = rcov, nObsEq = n,
         solvetol = solvetol )
      yCov <- .calcRCov( Y, rcovformula = 0, nObsEq = n, centered = TRUE,
         solvetol = solvetol )
      residCovInv <- solve( rcov, tol = solvetol )
      denominator <- 0
      for( i in 1:G ) {
         for( j in 1:G ) {
            denominator <- denominator + residCovInv[ i, j ] * yCov[ i, j ] * n[1]
         }
      }
      mcelr2 <- 1 - ( rtOmega %*% resids ) / denominator
#       # second formula from Greene (2003, p. 345)
#        yCov <- sum(diag(.calcRCov( Y, rcovformula = 0, nObsEq = n, centered = TRUE,
#           solvetol = solvetol )))
#        yCov <- drop( t(Y-mean(Y)) %*% (Y-mean(Y)) / sum(n) )
#       yCov <- .calcRCov( Y, rcovformula = 1, nObsEq = n, nCoefEq=rep(1,G), 
#          centered = TRUE, solvetol = solvetol )
#        mcelr2 <- 1 - G / sum( diag( solve( rcov, tol = solvetol ) * yCov ) )
  } else {
     mcelr2 <- NA
  }

  b              <- c(b)
  names(b)       <- xnames
  se             <- c(se)
  names(se)      <- xnames
  t              <- c(t)
  names(t)       <- xnames
  prob           <- c(prob)
  names(prob)    <- xnames
  colnames( bcov ) <- xnames
  rownames( bcov ) <- xnames
  colnames( rcov ) <- eqnlabels
  rownames( rcov ) <- eqnlabels

  ## build the "return" structure for the whole system
  results$method  <- method
  results$g       <- G              # number of equations
  results$n       <- N              # total number of observations
  results$k       <- K              # total number of coefficients
  results$ki      <- Ki             # total number of linear independent coefficients
  results$df      <- N - Ki         # dewgrees of freedom of the whole system
  results$b       <- b              # all estimated coefficients
  results$bt      <- bt             # transformed vector of estimated coefficients
  results$se      <- se             # standard errors of estimated coefficients
  results$t       <- t              # t-values of estimated coefficients
  results$p       <- prob           # p-values of estimated coefficients
  results$bcov    <- bcov           # coefficients covariance matrix
  results$btcov   <- btcov          # covariance matrix for transformed coeff. vector
  results$rcov    <- rcov           # residual covarance matrix
  results$drcov   <- drcov          # determinant of residual covarance matrix
  results$olsr2   <- olsr2          # R-squared value of the equation system
  results$iter    <- iter           # residual correlation matrix
  results$y       <- y              # vector of all (stacked) endogenous variables
  results$x       <- X              # matrix of all (diagonally stacked) regressors
  results$resids  <- resids         # vector of all (stacked) residuals
  results$data    <- alldata        # data frame for all data used in the system
  if( method == "2SLS" | method == "W2SLS" | method == "3SLS" |
      method == "W3SLS" ) {
    results$h       <- H            # matrix of all (diagonally stacked) instr. variables
    results$xHat    <- Xf           # matrix of "fitted" regressors
  }
  if( method == "SUR" | method == "WSUR" | method == "3SLS" |
      method == "W3SLS" ) {
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
  class(results)  <- "systemfit"

  results
}


## print the (summary) results that belong to the whole system
summary.systemfit <- function(object,...) {
   object$coef <- cbind( object$b, object$se, object$t, object$p )
   colnames( object$coef ) <- c( "Estimate", "Std. Error",
      "t value", "Pr(>|t|)" )
   object$rcor <- cor( residuals( object ) )
   dimnames( object$rcor ) <- dimnames( object$rcov )
   class( object ) <- "summary.systemfit"
   return( object )
}

## print summary results of the whole system
print.summary.systemfit <- function( x,
      digits = max( 3, getOption("digits") - 1 ), ... ) {

  table <- NULL
  labels <- NULL

  cat("\n")
  cat("systemfit results \n")
  cat("method: ")
  if(!is.null(x$iter)) if(x$iter>1) cat("iterated ")
  cat( paste( x$method, "\n\n"))
  if(!is.null(x$iter)) {
    if(x$iter>1) {
      if(x$iter<x$maxiter) {
        cat( paste( "convergence achieved after",x$iter,"iterations\n\n" ) )
      } else {
        cat( paste( "warning: convergence not achieved after", x$iter,
                    "iterations\n\n" ) )
      }
    }
  }
  for(i in 1:x$g) {
    row <- NULL
    row <- cbind( round( x$eq[[i]]$n,     digits ),
                  round( x$eq[[i]]$df,    digits ),
                  round( x$eq[[i]]$ssr,   digits ),
                  round( x$eq[[i]]$mse,   digits ),
                  round( x$eq[[i]]$rmse,  digits ),
                  round( x$eq[[i]]$r2,    digits ),
                  round( x$eq[[i]]$adjr2, digits ))
    table  <- rbind( table, row )
    if( is.null( x$eq[[i]]$eqnlabel ) ) {
      labels <- rbind( labels, paste( "equation", i ) )
    } else {
      labels <- rbind( labels, x$eq[[i]]$eqnlabel )
    }
  }
  rownames(table) <- c( labels )
  colnames(table) <- c("N","DF", "SSR", "MSE", "RMSE", "R2", "Adj R2" )

  ##print.matrix(table, quote = FALSE, right = TRUE )
  ##prmatrix(table, quote = FALSE, right = TRUE )
  print(table, quote = FALSE, right = TRUE, digits = digits )

  cat("\n")

  if(!is.null(x$rcovest)) {
    cat("The covariance matrix of the residuals used for estimation\n")
    rcov <- x$rcovest
    rownames(rcov) <- labels
    colnames(rcov) <- labels
    print( rcov, digits = digits )
    cat("\n")
    if( min(eigen( x$rcov, only.values=TRUE)$values) < 0 ) {
      cat("warning: this covariance matrix is NOT positive semidefinit!\n")
      cat("\n")
    }
  }

  cat("The covariance matrix of the residuals\n")
  rcov <- x$rcov
  rownames(rcov) <- labels
  colnames(rcov) <- labels
  print( rcov, digits = digits )
  cat("\n")

  cat("The correlations of the residuals\n")
  rcor <- x$rcor
  rownames(rcor) <- labels
  colnames(rcor) <- labels
  print( rcor, digits = digits )
  cat("\n")

  cat("The determinant of the residual covariance matrix: ")
  cat( formatC( x$drcov, digits = digits, width = -1 ) )
  cat("\n")

  cat("OLS R-squared value of the system: ")
  cat( formatC( x$olsr2, digits = digits, width = -1 ) )
  cat("\n")

  if(!is.null(x$mcelr2)) {
    cat("McElroy's R-squared value for the system: ")
    cat( formatC( x$mcelr2, digits = digits, width = -1 ) )
    cat("\n")
  }
  ## now print the individual equations
  for(i in 1:x$g) {
      print( summary( x$eq[[i]] ), digits = digits )
  }
  invisible( x )
}

## print a few results of the whole system
print.systemfit <- function( x,
      digits = max( 3, getOption("digits") - 1 ),... ) {

   cat("\n")
   cat("systemfit results \n")
   cat("method: ")
   if(!is.null(x$iter)) if(x$iter>1) cat("iterated ")
   cat( paste( x$method, "\n\n"))
   if(!is.null(x$iter)) {
      if(x$iter>1) {
         if(x$iter<x$maxiter) {
            cat( paste( "convergence achieved after",x$iter,"iterations\n\n" ) )
         } else {
            cat( paste( "warning: convergence not achieved after", x$iter,
                        "iterations\n\n" ) )
         }
      }
   }
   cat( "Coefficients:\n" )
   print( x$b, digits = digits )
   invisible( x )
}

## print the (summary) results for a single equation
summary.systemfit.equation <- function(object,...) {
   object$coef <- cbind( object$b, object$se, object$t, object$p )
   colnames( object$coef ) <- c( "Estimate", "Std. Error",
      "t value", "Pr(>|t|)" )
   class( object ) <- "summary.systemfit.equation"
   return( object )
}


## print summary results for a single equation
print.summary.systemfit.equation <- function( x,
   digits = max( 3, getOption("digits") - 1 ), ... ) {

  cat("\n")
  if( is.null( x$eqnlabel ) ) {
    cat( x$method, " estimates for equation ", x$i, "\n", sep = "" )
  } else {
    cat( x$method, " estimates for '", x$eqnlabel,
         "' (equation ", x$i, ")\n", sep = "" )
  }

  cat("Model Formula: ")
  print(x$formula)
  if(!is.null(x$inst)) {
    cat("Instruments: ")
    print(x$inst)
  }
  cat("\n")

  Signif <- symnum(x$p, corr = FALSE, na = FALSE,
                   cutpoints = c( 0, 0.001, 0.01, 0.05, 0.1, 1 ),
                   symbols   = c( "***", "**", "*", "." ," " ))

  table <- cbind(round( x$b,  digits ),
                 round( x$se, digits ),
                 round( x$t,  digits ),
                 round( x$p,  digits ),
                 Signif)

  rownames(table) <- names(x$b)
  colnames(table) <- c("Estimate","Std. Error","t value","Pr(>|t|)","")

  ##print.matrix(table, quote = FALSE, right = TRUE )
  ##prmatrix(table, quote = FALSE, right = TRUE )
  print(table, quote = FALSE, right = TRUE, digits = digits )
  cat("---\nSignif. codes: ",attr(Signif,"legend"),"\n")

  cat(paste("\nResidual standard error:", round(x$s, digits),
            "on", x$df, "degrees of freedom\n"))
            # s ist the variance, isn't it???

  cat( paste( "Number of observations:", round(x$n, digits),
              "Degrees of Freedom:", round(x$df, digits),"\n" ) )

  cat( paste( "SSR:",     round(x$ssr,    digits),
              "MSE:", round(x$mse, digits),
              "Root MSE:",   round(x$rmse,  digits), "\n" ) )

  cat( paste( "Multiple R-Squared:", round(x$r2,    digits),
              "Adjusted R-Squared:", round(x$adjr2, digits),
              "\n" ) )
  cat("\n")
  invisible( x )
}

## print a few results for a single equation
print.systemfit.equation <- function( x,
      digits = max( 3, getOption("digits") - 1 ), ... ) {

   cat("\n")
   if( is.null( x$eqnlabel ) ) {
      cat( x$method, " estimates for equation ", x$i, "\n", sep = "" )
   } else {
      cat( x$method, " estimates for '", x$eqnlabel,
            "' (equation ", x$i, ")\n", sep = "" )
   }

   cat("Model Formula: ")
   print(x$formula)
   if(!is.null(x$inst)) {
      cat("Instruments: ")
      print(x$inst)
   }

   cat("\nCoefficients:")
   print( x$b, digits = digits )
   invisible( x )
}


## calculate predicted values, its standard errors and the prediction intervals
predict.systemfit <- function( object, data=object$data,
                               se.fit=FALSE, se.pred=FALSE,
                               interval="none", level=0.95, ... ) {
   attach(data); on.exit( detach( data ) )

   predicted <- data.frame( obs=seq( nrow( data ) ) )
   colnames( predicted ) <- as.character( 1:ncol( predicted ) )
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
   if( object$method == "SUR" | object$method == "WSUR" |
       object$method == "3SLS" | object$method == "W3SLS") {
      if( se.fit | interval == "confidence" ) {
         ycovc <- X %*% object$bcov %*% t(X)
      }
      if( se.pred | interval == "prediction" ) {
         ycovp <- X %*% object$bcov %*% t(X) + object$rcov %x% diag(1,n[1],n[1])
      }
   }
   for(i in 1:g) {
      # fitted values
      Yi <- Y[(1+sum(n[1:i])-n[i]):sum(n[1:i]),]
      predicted <- cbind( predicted, Yi )
      names( predicted )[ length( predicted ) ] <- paste( "eq", as.character(i),
                                                          ".pred", sep="" )
      # calculate variance covariance matrices
      if( se.fit | interval == "confidence" ) {
         if( object$method == "SUR" | object$method == "WSUR" |
             object$method == "3SLS" | object$method== "W3SLS" ) {
            ycovci <- ycovc[ ( 1 + sum( n[1:i] ) - n[i] ) : sum( n[1:i] ),
                             ( 1 + sum( n[1:i] ) - n[i] ) : sum( n[1:i] ) ]
         } else {
            ycovci <- x[[i]] %*% object$eq[[i]]$bcov %*% t(x[[i]])
         }
      }
      if( se.pred | interval == "prediction" ) {
         if( object$method == "SUR" | object$method == "WSUR" |
             object$method == "3SLS" | object$method == "W3SLS" ) {
            ycovpi <- ycovp[ ( 1 + sum( n[1:i] ) - n[i] ) : sum( n[1:i] ),
                            ( 1 + sum( n[1:i] ) - n[i] ) : sum( n[1:i] ) ]
         } else {
            ycovpi <- x[[i]] %*% object$eq[[i]]$bcov %*% t(x[[i]]) +
                                 object$eq[[i]]$s2
         }
      }
      # standard errors of fitted values
      if( se.fit ) {
         if(nrow(data)==1) {
            predicted <- cbind( predicted, sqrt( ycovci ) )
         } else {
            predicted <- cbind( predicted, sqrt( diag( ycovci ) ) )
         }
         names( predicted )[ length( predicted ) ] <-
            paste( "eq", as.character(i), ".se.fit", sep="" )
      }
      # standard errors of prediction
      if( se.pred ) {
         if(nrow(data)==1) {
            predicted <- cbind( predicted, sqrt( ycovpi ) )
         } else {
            predicted <- cbind( predicted, sqrt( diag( ycovpi ) ) )
         }
         names( predicted )[ length( predicted ) ] <-
            paste( "eq", as.character(i), ".se.pred", sep="" )
      }

      # confidence intervals
      if( interval == "confidence" ) {
         if( object$probdfsys ) {
            tval   <- qt( 1 - ( 1- level )/2, object$df )
         } else {
            tval   <- qt( 1 - ( 1- level )/2, object$eq[[i]]$df )
         }
         if( nrow(data)==1 ) {
            seci    <- sqrt( ycovci )
         } else {
            seci    <- sqrt( diag( ycovci ) )
         }
         predicted <- cbind( predicted, Yi - ( tval * seci ) )
         names( predicted )[ length( predicted ) ] <-
            paste( "eq", as.character(i), ".lwr", sep="" )
         predicted <- cbind( predicted, Yi + ( tval * seci ) )
         names( predicted )[ length( predicted ) ] <-
            paste( "eq", as.character(i), ".upr", sep="" )
      }
      # prediction intervals
      if( interval == "prediction" ) {
         if( object$probdfsys ) {
            tval   <- qt( 1 - ( 1- level )/2, object$df )
         } else {
            tval   <- qt( 1 - ( 1- level )/2, object$eq[[i]]$df )
         }
         if(nrow(data)==1) {
            sepi <- sqrt( ycovpi )
         } else {
            sepi <- sqrt( diag( ycovpi ) )
         }
         predicted <- cbind( predicted, Yi - ( tval * sepi ) )
         names( predicted )[ length( predicted ) ] <-
            paste( "eq", as.character(i), ".lwr", sep="" )
         predicted <- cbind( predicted, Yi + ( tval * sepi ) )
         names( predicted )[ length( predicted ) ] <-
            paste( "eq", as.character(i), ".upr", sep="" )
      }
   }
   predicted[ 2: length( predicted ) ]
}

## calculate predicted values, its standard errors and the prediction intervals
predict.systemfit.equation <- function( object, data=object$data, ... ) {
   attach( data ); on.exit( detach( data ) )
   x <-  model.matrix( object$formula )
   predicted <- drop( x %*% object$b )
   predicted
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


## return all coefficients
coef.systemfit <- function( object, ... ) {
   object$b
}

## return all coefficients, std.errors, t-values and p-values
coef.summary.systemfit <- function( object, ... ) {
   object$coef
}

## return the coefficients of a single equation
coef.systemfit.equation <- function( object, ... ) {
   object$b
}

## return coefficients, std.errors, t-values and p-values of a single equation
coef.summary.systemfit.equation <- function( object, ... ) {
   object$coef
}

## return all residuals
residuals.systemfit <- function( object, ... ) {
   if( is.null( colnames( object$rcov ) ) ) {
      eqNames <- paste( "eq", c( 1:object$g ) )
   } else {
      eqNames <- colnames( object$rcov )
   }
   residuals <- data.frame( obsNo = c( 1:length( object$eq[[1]]$residuals ) ) )
   for( i in 1:object$g ) {
      residuals[[ eqNames[ i ] ]] <- object$eq[[i]]$residuals
   }
   residuals$obsNo <- NULL
   return( residuals )
}

## return residuals of a single equation
residuals.systemfit.equation <- function( object, ... ) {
   object$residuals
}

## return the variance covariance matrix of the coefficients
vcov.systemfit <- function( object, ... ) {
   object$bcov
}

## return the variance covariance matrix of the coefficients of a single equation
vcov.systemfit.equation <- function( object, ... ) {
   object$bcov
}

## return the variance covariance matrix of the coefficients
confint.systemfit <- function( object, parm = NULL, level = 0.95, ... ) {
   a <- ( 1 - level ) / 2
   a <- c( a, 1 - a )
   pct <- paste( round( 100 * a, 1 ), "%" )
   ci <- array( NA, dim = c( length( object$b ), 2),
            dimnames = list( names( object$b ), pct ) )
   j <- 1
   for( i in 1:object$g ) {
      object$eq[[i]]$dfSys <- object$df
      object$eq[[i]]$probdfsys <- object$probdfsys
      ci[ j:(j+object$eq[[ i ]]$k-1), ] <- confint( object$eq[[ i ]] )
      j <- j + object$eq[[ i ]]$k
   }
   class( ci ) <- "confint.systemfit"
   ci
}

## return the variance covariance matrix of the coefficients of a single equation
confint.systemfit.equation <- function( object, parm = NULL, level = 0.95, ... ) {
   a <- ( 1 - level ) / 2
   a <- c( a, 1 - a )
   pct <- paste( round( 100 * a, 1 ), "%" )
   ci <- array( NA, dim = c( length( object$b ), 2),
            dimnames = list( names( object$b ), pct ) )
   if( object$ probdfsys ) {
      fac <- qt( a, object$dfSys )
   } else {
      fac <- qt( a, object$df )
   }
   ci[] <- object$b + object$se %o% fac
   class( ci ) <- "confint.systemfit"
   ci
}

## print the confidence intervals of the coefficients
print.confint.systemfit <- function( x, digits = 3, ... ) {
   print( unclass( round( x, digits = digits, ...) ) )
   invisible(x)
}

## return the fitted values
fitted.systemfit <- function( object, ... ) {
   fitted <- array( NA, c( length( object$eq[[1]]$fitted ), object$g ) )
   colnames( fitted ) <- as.character( 1:ncol( fitted ) )
   for(i in 1:object$g )  {
      fitted[ , i ]           <- object$eq[[ i ]]$fitted
      colnames( fitted )[ i ] <- paste( "eq", as.character(i), sep="" )
   }
   fitted
}

## return the fitted values of e single euation
fitted.systemfit.equation <- function( object, ... ) {
   object$fitted
}




