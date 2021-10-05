cvamMlogit.fit <- function( x, y, baseline=1L, criterion=1e-06,
   iterMax=25L ) {
   stopifnot( NROW(x) == NROW(y) )
   stopifnot( NCOL(y) >= 2L )
   n <- NROW(x)
   p <-  NCOL(x)
   r <- NCOL(y)
   storage.mode(y) <- "double"
   storage.mode(x) <- "double"
   baseline <- as.integer(baseline)[1L]
   criterion <- as.double(criterion)[1L]
   iterMax <- as.integer(iterMax)[1L]
   #--------------------------------------------
   # create a matrix for holding message codes
   msg.len.max <- 40L
   msg.codes <- matrix( 0L, msg.len.max, 17L )
   #--------------------------------------------
   tmp <- .Fortran("cvam_mlogit",
      n = n,
      p = p,
      r = r,
      x = x,
      y = y,
      baseline = baseline,
      iterMax = iterMax,
      criterion = criterion,
      iter = integer(1L),
      convergedInt = integer(1L),
      loglik = numeric(1L),
      score = numeric( p*(r-1L) ),
      hess = matrix( numeric(1L), p*(r-1L), p*(r-1L) ),
      beta = matrix( numeric(1L), p, r ),
      betaVec = numeric( p*(r-1L) ),
      vhatBetaVec = matrix( numeric(1L), p*(r-1L), p*(r-1L) ),
      piMat = matrix( numeric(1L), n, r ),
      # messaging
      status = integer(1L),
      msg.len.max = msg.len.max,
      msg.codes = msg.codes,
      msg.len.actual = integer(1L),
      PACKAGE = "cvam" )
      #--------------------------------------------
      # display message from Fortran, if present
      msg.lines <- .msgCvam( tmp$msg.codes, tmp$msg.len.actual )
      if( is.null( msg.lines ) ){
         msg <- "OK"
      } else{
         msg <- paste0( msg.lines, collapse="\n" )
      }
      msg <- paste( msg, "\n", sep="")
      if( msg!= "OK\n" ) cat( paste("Note: ", msg, sep="") )
   #--------------------------------------------
   if( tmp$status != 0 ) stop( gettext( 
      "Procedure aborted" ), domain = NA )
   #--------------------------------------------
   tmp$converged <- as.logical(tmp$convergedInt)
   tmp
}

cvamMlogit <- function( form, data, baseline=1L,
   prior = c("none", "DAP"), priorFreq = NULL, 
   criterion=1e-06, iterMax=25L ) {
   #--------------------------------------------
   # create model frame, extract x and y
   mc <- match.call( expand.dots=FALSE )
   mc[[1]] <- quote(stats::model.frame)
   mc$baseline <- mc$prior <- mc$priorFreqTot <-
      mc$criterion <- mc$iterMax <- NULL
   m <- match( c("form", "data"), names(mc), nomatch=0L )
   mc <- mc[ c(1L,m) ]
   names(mc)[2] <- "formula"
   mc$na.action <- as.name("na.fail")
   mc$drop.unused.levels <- FALSE
   mF <- eval( mc, parent.frame() )
   y <- model.response(mF)
   stopifnot( NCOL(y) >= 2L )
   x <- model.matrix(form, mF)
   #--------------------------------------------
   prior <- match.arg( prior )
   if( prior == "DAP" ) {
      if( is.null(priorFreq) ) {
         priorFreqTot <- as.double( NCOL(x) * ( NCOL(y) - 1L ) )
         interceptOnlyFit <- cvamMlogit.fit( rep(1,NROW(x)),
            y, baseline, criterion, iterMax )
         if( ! interceptOnlyFit$converged ) stop( gettext(
            "Intercept-only model failed to converge" ), domain = NA )
         alpha <- interceptOnlyFit$piMat[1L,]
         stopifnot( is.null(mF$`_freq`) )
         mFGrouped <- mF
         mFGrouped$`_freq` <- 1
         form2 <- form
         form2[[2L]] <- as.symbol("_freq")
         mFGrouped <- aggregate(form2, data=mFGrouped, FUN=sum )
         nCovPatt <- NROW(mFGrouped)
         priorFreqPerPatt <- priorFreqTot / nCovPatt
         stopifnot( is.null(mF$`_ROW`) )
         mF$`_ROW` <- 1:NROW(mF)
         mF <- merge( mF, mFGrouped, all.x=TRUE, all.y=FALSE, sort=FALSE )
         mF <- mF[ order(mF$`_ROW`), ]
         freq <- mF$`_freq`
         priorFreq <- matrix( alpha, NROW(y), NCOL(y), byrow=TRUE ) *
             ( priorFreqPerPatt / freq )
      } else {
         stopifnot( NROW(priorFreq) == NROW(y) )
         stopifnot( NCOL(priorFreq) == NCOL(y) )
         stopifnot( all( priorFreq >=0 ) )
      }
   } else {
      priorFreq <- matrix(0, NROW(y), NCOL(y) )
   }
   #--------------------------------------------
   result <- cvamMlogit.fit( x, y + priorFreq, baseline, criterion, iterMax )
   result$prior <- prior
   result$priorFreq <- if( prior == "DAP" ) priorFreq else NULL
   result$alpha <- if( prior == "DAP" ) alpha else NULL
   result$y <- if( prior == "DAP" ) result$y - priorFreq else result$y
   result$fitted <- result$piMat * apply(result$y, 1, sum )
   result$fittedWithPrior <- if( prior == "DAP" )
      result$piMat * apply(result$y+result$priorFreq, 1, sum ) else NULL
   result$logP <- result$loglik
   if( prior == "DAP" ) {
      tmp <- result$priorFreq * log( result$piMat )
      tmp[ result$priorFreq==0 ] <- 0
      logPrior <- sum(tmp)
      result$loglik <- result$loglik - logPrior
   }
   return(result)
}

.piMat <- function( xBeta, likMat=NULL ) {
   #--- convert xBeta to probs
   shift <- apply(xBeta, 1, max)
   piMat <- exp(xBeta - shift)
   piMat <- piMat / apply(piMat, 1, sum )
   #--- multiply by likMat and normalize
   if( !is.null(likMat) ) piMat <- piMat * likMat
   piMat <- piMat / apply(piMat, 1, sum )
   piMat
}

fitLCPrev <- function( form, data, likMat, freq,
   baseline=NULL, prior=c("none", "DAP"), priorFreq=NULL,
   criterion=1e-06, iterMaxEM=500L, iterMaxNR=25L ) {
   #-------------------------
   # create model frame, extract model.matrix and model.response
   mc <- match.call( expand.dots=FALSE )
   mc[[1]] <- quote(stats::model.frame)
   m <- match( c("form", "data", "freq"), names(mc), nomatch=0L )
   mc <- mc[ c(1L,m) ]
   names(mc)[2] <- "formula"
   mc$na.action <- as.name("na.fail")
   mc$drop.unused.levels <- FALSE
   mF <- eval( mc, parent.frame() )
   x <- model.matrix(form, mF)
   L <- model.response(mF)
   if( ! is.latentFactor(L) ) stop( gettext(
      "Response variable in formula is not a latent factor" ), domain = NA )
   if( is.null( mF$`(freq)` ) ) {
      freq <- rep(1L, NROW(mF))
      freqSupplied <- FALSE
   } else {
      freq <- mF$`(freq)`
      freqSupplied <- TRUE
   }
   storage.mode(freq) <- "double"
   stopifnot( all(freq >= 0) )
   #-------------------------
   # check likMat and baseline
   if( NROW(likMat) != NROW(x) ) stop( gettextf(
      "likMat has incorrect number of rows, should be %i", NROW(x)),
      domain = NA )
   if( NCOL(likMat) != nBaseLevels(L) ) stop( gettextf(
      "likMat has incorrect number of columns, should be %i", nBaseLevels(L)),
      domain = NA )
   if( ! setequal( colnames(likMat), baseLevels(L) ) ) stop( gettext(
      "colnames(likMat) are not the baseLevels of the latent factor" ),
       domain = NA )
   if( is.null(baseline) ) {
      baseline <- baseLevels(L)[[1L]]
   } else {
      baseline <- as.character(baseline)[[1L]]
      if( ! ( baseline %in% baseLevels(L) ) ) stop( gettext(
         "baseline is not one of the baseLevels of the latent factor" ),
         domain = NA )
   }
   m <- match( baseLevels(L), colnames(likMat) )
   likMat <- data.matrix( likMat[,m] )
   baseline <- match( baseline, baseLevels(L) )
   #-------------------------
   # handle prior
   prior <- match.arg( prior )
   if( prior == "DAP" ) {
      if( is.null(priorFreq) ) {
         priorFreqTot <- as.double( NCOL(x) * ( nBaseLevels(L) - 1L ) )
         xInt <- rep(1,NROW(x))
         betaNew <- matrix(0, 1L, nBaseLevels(L) )
         converged <- FALSE
         iter <- 0L
         while( ( ! converged ) & ( iter < iterMaxEM ) ) {
	    iter <- iter + 1L
            beta <- betaNew
            Lhat <- .piMat( xInt %*% beta, likMat ) * freq
            fit <- cvamMlogit.fit( xInt, Lhat, baseline, criterion, iterMaxNR )
            if( ! fit$converged ) stop( gettext(
               "Intercept-only model failed to converge" ), domain = NA )
            betaNew[] <- fit$beta
            converged <- all( abs(betaNew-beta) <= criterion )
         }
         if( ! converged ) stop( gettext(
            "Intercept-only model failed to converge" ), domain = NA )
         alpha <- fit$piMat[1L,]
         stopifnot( is.null(mF$`_freq`) )
         mFGrouped <- mF
         mFGrouped$`_freq` <- 1
         form2 <- form
         form2[[2L]] <- as.symbol("_freq")
         mFGrouped <- aggregate(form2, data=mFGrouped, FUN=sum )
         nCovPatt <- NROW(mFGrouped)
         priorFreqPerPatt <- priorFreqTot / nCovPatt
         stopifnot( is.null(mF$`_ROW`) )
         mF$`_ROW` <- 1:NROW(mF)
         mF <- merge( mF, mFGrouped, all.x=TRUE, all.y=FALSE, sort=FALSE )
         mF <- mF[ order(mF$`_ROW`), ]
         pfreq <- mF$`_freq`
         priorFreq <- matrix( alpha, NROW(x), nBaseLevels(L), byrow=TRUE ) *
             ( priorFreqPerPatt / pfreq )
      } else {
         stopifnot( NROW(priorFreq) == NROW(x) )
         stopifnot( NCOL(priorFreq) == nBaseLevels(L) )
         stopifnot( all( priorFreq >=0 ) )
      }
   } else {
      priorFreq <- matrix(0, NROW(x), nBaseLevels(L) )
   }
   #-------------------------
   betaNew <- matrix(0, NCOL(x), nBaseLevels(L) )
   rownames(betaNew) <- colnames(x)
   colnames(betaNew) <- baseLevels(L)
   converged <- FALSE
   iter <- 0L
   aborted <- FALSE
   while( ( ! converged ) & ( iter < iterMaxEM ) ) {
      iter <- iter + 1L
      beta <- betaNew
      Lhat <- .piMat( x %*% beta, likMat ) * freq
      fit <- cvamMlogit.fit( x, Lhat+priorFreq, baseline,
         criterion, iterMaxNR )
      if( ! fit$converged ) {
         aborted <- TRUE
         break
      }
      betaNew[] <- fit$beta
      converged <- all( abs(betaNew-beta) <= criterion )
   }
   #-------------------------
   beta <- betaNew
   betaVec <- fit$betaVec
   pi <- .piMat( x %*% beta )
   piStar <- .piMat( x %*% beta, likMat )
   Lhat <- piStar * freq
   #-------------------------
   loglik.derivs <- .LCPrev.loglik.derivs( x, betaVec, likMat, freq, 
       baseline)
   logP <- loglik <- loglik.derivs$loglik
   score <- loglik.derivs$score
   hess <-  loglik.derivs$hess
   if( prior == "DAP" ) {
      logprior.derivs <-
         .Mlogit.loglik.derivs( x, priorFreq, betaVec, baseline )
      logP <- logP + logprior.derivs$loglik
      score <- score + logprior.derivs$score
      hess <-  hess + logprior.derivs$hess
   }
   tmp <- try( chol(-hess), silent=TRUE )
   if( inherits(tmp, "try-error") ) {
      message( gettext(
         "logP not concave, standard errors not available" ), domain = NA )
      vhatBetaVec <- NULL
   } else vhatBetaVec <- solve(-hess)
   #-------------------------
   result <- list(
      beta = betaNew,
      pi = pi,
      piStar = piStar,
      Lhat = Lhat,
      prior = prior,
      priorFreq = if( prior == "DAP" ) priorFreq else NULL,
      alpha = if( prior == "DAP" ) alpha else NULL,
      x = x,
      freq = freq,
      betaVec = betaVec,
      baseline = baseline,
      likMat = likMat,
      aborted = aborted,
      iter = iter,
      converged = converged,
      loglik = loglik,
      logP = logP,
      score = score,
      hess = hess,
      vhatBetaVec = vhatBetaVec )
   return(result)
}

.Mlogit.loglik.derivs <- function( x, y, betaVec, baseline=1L ) {
   stopifnot( NROW(x) == NROW(y) )
   stopifnot( NCOL(y) >= 2L )
   n <- NROW(x)
   p <-  NCOL(x)
   r <- NCOL(y)
   stopifnot( length(betaVec) == p*(r-1L) )
   storage.mode(betaVec) <- storage.mode(y) <- storage.mode(x) <- "double"
   baseline <- as.integer(baseline)[1L]
   #--------------------------------------------
   # create a matrix for holding message codes
   msg.len.max <- 40L
   msg.codes <- matrix( 0L, msg.len.max, 17L )
   #--------------------------------------------
   tmp <- .Fortran("cvam_mlogit_loglik_derivs",
      n = n,
      p = p,
      r = r,
      x = x,
      y = y,
      baseline = baseline,
      betaVec <- betaVec,
      loglik = numeric(1L),
      score = numeric( p*(r-1L) ),
      hess = matrix( numeric(1L), p*(r-1L), p*(r-1L) ),
      # messaging
      status = integer(1L),
      msg.len.max = msg.len.max,
      msg.codes = msg.codes,
      msg.len.actual = integer(1L),
      PACKAGE = "cvam" )
      #--------------------------------------------
      # display message from Fortran, if present
      msg.lines <- .msgCvam( tmp$msg.codes, tmp$msg.len.actual )
      if( is.null( msg.lines ) ){
         msg <- "OK"
      } else{
         msg <- paste0( msg.lines, collapse="\n" )
      }
      msg <- paste( msg, "\n", sep="")
      if( msg!= "OK\n" ) cat( paste("Note: ", msg, sep="") )
   #--------------------------------------------
   if( tmp$status != 0 ) stop( gettext( 
      "Procedure aborted" ), domain = NA )
   #--------------------------------------------
   list(
      loglik = tmp$loglik,
      score = tmp$score,
      hess = tmp$hess )
}

.LCPrev.loglik.derivs <- function( x, betaVec, likMat, freq=NULL, 
      baseline=1L ) {
   stopifnot( NROW(x) == NROW(likMat) )
   stopifnot( NCOL(likMat) >= 2L )
   n <- NROW(x)
   p <-  NCOL(x)
   r <- NCOL(likMat)
   stopifnot( length(betaVec) == p*(r-1L) )
   if( is.null(freq) ) freq <- rep(1, n)
   stopifnot( length(freq) == n )
   storage.mode(freq) <- storage.mode(betaVec) <-
      storage.mode(likMat) <- storage.mode(x) <- "double"
   baseline <- as.integer(baseline)[1L]
   #--------------------------------------------
   # create a matrix for holding message codes
   msg.len.max <- 40L
   msg.codes <- matrix( 0L, msg.len.max, 17L )
   #--------------------------------------------
   tmp <- .Fortran("cvam_lcprev_loglik_derivs",
      n = n,
      p = p,
      r = r,
      x = x,
      likMat = likMat,
      freq = freq,
      baseline = baseline,
      betaVec <- betaVec,
      loglik = numeric(1L),
      score = numeric( p*(r-1L) ),
      hess = matrix( numeric(1L), p*(r-1L), p*(r-1L) ),
      # messaging
      status = integer(1L),
      msg.len.max = msg.len.max,
      msg.codes = msg.codes,
      msg.len.actual = integer(1L),
      PACKAGE = "cvam" )
      #--------------------------------------------
      # display message from Fortran, if present
      msg.lines <- .msgCvam( tmp$msg.codes, tmp$msg.len.actual )
      if( is.null( msg.lines ) ){
         msg <- "OK"
      } else{
         msg <- paste0( msg.lines, collapse="\n" )
      }
      msg <- paste( msg, "\n", sep="")
      if( msg!= "OK\n" ) cat( paste("Note: ", msg, sep="") )
   #--------------------------------------------
   if( tmp$status != 0 ) stop( gettext( 
      "Procedure aborted" ), domain = NA )
   #--------------------------------------------
   list(
      loglik = tmp$loglik,
      score = tmp$score,
      hess = tmp$hess )
}
