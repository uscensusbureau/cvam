.isOneSided <- function( obj ) { inherits(obj, "formula") &&
   length(obj) == 2L }

cvam <- function( obj, ...) {
   # S3 generic function
   UseMethod("cvam")
}

cvam.default <- function( obj, ...) {
   stop( gettext(
      'First argument must be an object of class "formula" or "cvam"'),
      domain = NA )
   }

cvam.cvam <- function( obj, ...) {
   stop( gettext(
      "Not implemented yet"),
      domain = NA )
   }

cvamControl <- function( iterMaxEM = 500L, iterMaxNR = 50L,
   iterApproxBayes = 1L, imputeApproxBayes = FALSE,
   iterMCMC = 5000L, burnMCMC = 0L, thinMCMC = 1L, imputeEvery = 0L,
   saveProbSeries = FALSE,
   typeMCMC = c("DA","RWM"), tuneDA = c(10,.8,.8),
   tuneRWM = c(1000,.1),
   stuckLimit = 25L,
   startValDefault = c("center", "uniform"), startValJitter = 0,
   critEM = 1e-06, critNR = 1e-06, critBoundary = 1e-08, ncolMaxMM = 1000L,
   excludeAllNA = TRUE, critModelCheck=1e-08, confidence=.95, 
   probRound = TRUE, probDigits = 4L ) {
   stopifnot( ( iterMaxEM <- as.integer(iterMaxEM)[1L] ) >= 0L )
   stopifnot( ( iterMaxNR <- as.integer(iterMaxNR)[1L] ) >= 0L )
   stopifnot( ( iterApproxBayes <- as.integer(iterApproxBayes)[1L] ) >= 0L )
   imputeApproxBayes <- as.logical( imputeApproxBayes )[1L]
   stopifnot( ( iterMCMC <- as.integer(iterMCMC)[1L] ) >= 0L )
   stopifnot( ( burnMCMC <- as.integer(burnMCMC)[1L] ) >= 0L )
   stopifnot( ( thinMCMC <- as.integer(thinMCMC)[1L] ) >= 1L )
   stopifnot( ( imputeEvery <- as.integer(imputeEvery)[1L] ) >= 0L )
   saveProbSeries <- as.logical( saveProbSeries )[1L]
   typeMCMC <- match.arg( typeMCMC )
   tuneDA <- as.numeric( tuneDA )[1:3]
   stopifnot( tuneDA[1L] > 0 )
   stopifnot( tuneDA[3L] > 0 )
   tuneRWM <- as.numeric(tuneRWM)[1:2]
   stopifnot( tuneRWM[1L] > 0 )
   stopifnot( tuneRWM[2L] > 0 )
   stopifnot( ( stuckLimit <- as.integer(stuckLimit)[1L] ) > 0 )
   startValDefault <- match.arg( startValDefault )
   stopifnot( ( startValJitter <- as.double(startValJitter)[1L] ) >= 0 )
   stopifnot( ( critEM <- as.double(critEM)[1L] ) > 0 )
   stopifnot( ( critNR <- as.double(critNR)[1L] ) > 0 )
   stopifnot( ( critBoundary <- as.double(critBoundary)[1L] ) > 0 )
   stopifnot( ( ncolMaxMM <- as.integer(ncolMaxMM)[1L] ) >= 0L )
   excludeAllNA <- as.logical(excludeAllNA)[1L]
   stopifnot( ( critModelCheck <-
      as.double(critModelCheck)[1L] ) >= 0 )
   stopifnot( ( confidence <- as.double(confidence)[1L] ) > 0 )
   stopifnot( confidence < 1 )
   probRound <- as.logical(probRound)[1L]
   probDigits <- as.integer(probDigits)[1L]
   list(
      iterMaxEM = iterMaxEM,
      iterMaxNR = iterMaxNR,
      iterApproxBayes = iterApproxBayes,
      imputeApproxBayes = imputeApproxBayes,
      iterMCMC = iterMCMC,
      burnMCMC = burnMCMC,
      thinMCMC = thinMCMC,
      imputeEvery = imputeEvery,
      saveProbSeries = saveProbSeries,
      typeMCMC = typeMCMC,
      tuneDA = tuneDA,
      tuneRWM = tuneRWM,
      stuckLimit = stuckLimit,
      startValDefault = startValDefault,
      startValJitter = startValJitter,
      critEM = critEM,
      critNR = critNR,
      critBoundary = critBoundary,
      ncolMaxMM = ncolMaxMM,
      excludeAllNA = excludeAllNA,
      critModelCheck = critModelCheck,
      confidence = confidence,
      probRound = probRound,
      probDigits = probDigits )
}

.cvamModel <- function( obj, saturated ) {
   # process the model formula into a .cvamModel object
   if( ! .isOneSided(obj) ) stop( gettext(
      "Argument 'obj' is not a one-sided formula" ), domain=NA )
   aV <- all.vars(obj)
   aN <- setdiff( all.names(obj), c("~", "(", "+", "|", ":", "*") )
   badN <- setdiff( aN, aV )
   if( length(badN) > 0L ) stop( gettextf(
      "In model formula 'obj': symbol '%s' not allowed", badN[1L] ),
      domain = NA )
   badV <- intersect( aV, c("freq", "prob", "weights", "offset", "likVal") )
   if( length(badV) > 0L ) stop( gettextf(
      "In model formula 'obj': variable cannot be named '%s'", badV[1L] ),
         domain = NA )
   Form <- Formula::Formula(obj)
   nParts <- length( attr(Form, "rhs") )
      if( nParts > 2L ) stop( gettext( 
     "In formula 'obj': multiple '|' detected"), domain = NA )
   leftVars <- all.vars( formula(Form, rhs=1L ) )
   rightVars <- if( nParts == 2L ) 
      all.vars( formula(Form, rhs=2L) )  else character()
   if( length( setdiff(leftVars, rightVars) ) == 0L ) stop( gettext(
      "In formula 'obj': no variables are being modeled"),
      domain = NA )
   if( length(rightVars) > 0L ) {
      if( any( ! (rightVars %in% leftVars) ) ) stop( gettextf(
         "In formula 'obj': variable '%s' not found before '|'",
         rightVars[!(rightVars %in% leftVars)] ), domain = NA )
      part2 <- formula(Form, rhs=2L) 
      aN <- setdiff( all.names(part2), c("~", "(", "+") )
      aV <- all.vars(part2)
      aN <- setdiff( aN, aV )
      if( length(aN) > 0L ) stop( gettextf( 
         "In formula 'obj': symbol '%s' not allowed after '|'", 
         aN[1L] ), domain = NA )
   }
   vars <- unique( c(leftVars, rightVars) )
   fixed <- vars %in% rightVars
   o <- order(fixed)
   vars <- vars[o]
   fixed <- fixed[o]
   formulaStr <- deparse(obj, width.cutoff=500L)
   formulaStr <- paste("~ ", substring(formulaStr, 2L), sep="" )
   saturated <- as.logical( saturated )[1L]
   structure( formula(Form, rhs=1L),
      vars = vars,
      fixed = fixed,
      formulaStr = formulaStr,
      saturated = saturated,
      class = c(".cvamModel", class(obj)) )
}

.priorNugget <- function( obj, prior) {
   if( is.null( names(obj) ) ) stop( gettext(
      "A component of prior has no 'names'"), domain = NA )
   nam <- names(obj)
   if( any( is.na( nam ) ) ) stop( gettext(
      "A component of prior has missing values in 'names'"), domain = NA )
   nam <- trimws(nam)
   if( any( nchar(nam) == 0L ) ) stop( gettext(
      "A component of prior has one or more 'names' that are blank"),
      domain = NA )
   if( is.null(obj$freq) ) stop( gettext(
      "A component of prior has no 'freq'"), domain = NA )
   freq <- as.double(obj$freq)[1L]
   if( freq < 0 ) stop( gettext(
      "A 'freq' in prior is negative"), domain = NA )
   obj$freq <- NULL
   if( !all( sapply(obj, is.character) ) ) stop( gettext(
      "A value in prior does not have mode 'character'"), domain = NA )
   if( !all( sapply(obj, length) == 1L ) ) stop( gettext(
      "A value in prior does not have length 1"), domain = NA )
   vars <- names(obj) # could be length zero
   values <- as.character( unlist(obj) ) # could be length zero
   if( length(vars) > 0L ) {
      w <- !is.na(values)
      vars <- vars[w]
      values <- values[w]
   }
   if( length(vars) > 0 ) {
      structure( list( vars = vars, values = values, freq = freq),
         class = c( "priorNugget", class(obj) ) )
   } else NULL
}

cvamPrior <- function( obj=list(), flatten=0, ridge=0, intensity=1 ) {
   if( !is.list(obj) ) stop( gettext( 
      "Argument 'obj' is not a list" ), domain = NA )
   m <- match( FALSE, unlist( lapply( obj, is.list ) ), nomatch=0L )
   if( m != 0L ) stop( gettextf( 
      "'obj'[[%i]] is not a list", m ), domain = NA )
   flatten <- as.double( flatten )[1L]
   intensity <- as.double( intensity )[1L]
   if( flatten < 0 ) stop( gettext( 
      "Argument 'flatten' cannot be negative" ), domain = NA )
   if( intensity < 0 ) stop( gettext( 
      "Argument 'intensity' cannot be negative" ), domain = NA )
   obj <- lapply( obj, .priorNugget )
   if( length(obj) > 0 ) {
      if( any(unlist( lapply(obj, is.null) ) ) ) warning( gettext(
         "Prior nuggets containing no information were dropped" ),
         domain = NA )
      obj <- obj[ ! unlist(lapply(obj, is.null))  ]
   }
   stopifnot( ( ridge <- as.numeric(ridge)[1L] ) >= 0 )
   structure( 
      list( nuggets=obj, flatten=flatten, ridge=ridge, intensity=intensity ),
      class = c( "cvamPrior", class(obj) ) ) 
}

summary.cvamPrior <- function( object, showNuggets=TRUE, ...) {
   stopifnot( inherits(object, "cvamPrior" ) )
   showNuggets <- as.logical(showNuggets)[1L]
   sL <- object$nuggets 
   allVars <- unique( unlist( lapply( sL, `[[`, "vars" ) ) )
   dF <- matrix( as.character(NA), length(sL), length(allVars),
      dimnames = list( names(sL), allVars ) )
   for( i in seq_along(sL) ) {
      m <- match( sL[[i]]$vars, allVars )
      dF[i,m] <- sL[[i]]$values
   }
   dF <- data.frame(dF)
   dF$freq <- unlist( lapply( sL, `[[`, "freq" ) )
   object$displayFrame <- if( length(dF) > 0 ) dF else NULL
   object$showNuggets <- showNuggets
   structure( object,
      class = c("summary.cvamPrior", class(object) ) )
}

print.cvamPrior <- function(x, ...) {
   stopifnot( inherits(x, "cvamPrior") )
   cat( gettext( 'Object of class "cvamPrior", use "summary"' ), sep="\n")
   invisible()
}

print.summary.cvamPrior <- function(x, 
      showNuggets=x$showNuggets, ...) {
   stopifnot( inherits(x, "summary.cvamPrior") )
   if( !is.null( x$displayFrame ) & showNuggets ) {
      cat( "Prior information nuggets:", sep="\n")
      print.data.frame( x$displayFrame, ...)
      cat("\n")
   }
   strA <- format( c("Flattening frequency =",
      "Total nuggets + flattening =",
      "Ridge factor =",
      "Intensity factor ="), justify="right")
   tF <- x$flatten + sum(x$displayFrame$freq)
   strB <- format( c( x$flatten, tF, x$ridge, x$intensity) )
   cat( paste(strA, strB), sep="\n")
   invisible()
}

.cvamPriorSummaries <- function( prior, pF ){
   pFreq <- unlist( lapply( prior$nuggets, `[[`, "freq" ) ) # NULL if no nuggets
   if( is.null(pFreq) ){
      pFreq <- numeric(0L)
      pFreqInt <- integer(0L)
   }
   for(i in seq_along(prior$nuggets) ) {
      pF[i,] <- NA
      vars <- prior$nuggets[[i]]$vars
      values <- prior$nuggets[[i]]$values
      for(j in seq_along(vars) ) {
         if( vars[j] %in% names(pF) ) {
            pFj <- match( vars[j], names(pF) )
            if( ! values[j] %in% levels( pF[[pFj]] ) ) stop( gettextf(
               "In 'prior': '%s' is not a level of '%s'",
               values[j], vars[j] ), domain = NA )
            pF[i, pFj] <- values[j]
         }
      }
   }
   if( length(prior$nuggets) > 0 ){
      allMiss <- apply( is.na(pF), 1, all )
      if( any(allMiss) ) {
         warning( gettextf(
            paste("Prior nuggets with no information about",
            "modeled variables were dropped", sep="\n")), domain=NA )
         pF <- pF[!allMiss,]
         pFreq <- pFreq[!allMiss]
         
      }
   }
   pFreq <- pFreq * prior$intensity
   pF$freq <- as.double( pFreq )
   pF$freqInt <- as.integer( pFreq )
   list( pF = pF, flatten = prior$flatten * prior$intensity,
      ridge = prior$ridge * prior$intensity ) 
}

cvam.formula <- function( obj, data, freq, prior=cvamPrior(), 
   method=c("EM", "MCMC", "approxBayes", "mfSeen",
      "mfTrue", "mfPrior", "modelMatrix"),
   control=list(), omitData=FALSE, saturated=FALSE, 
   modelMatrix = NULL, offset = NULL, strZero = NULL,
   startVal = NULL, estimate = NULL, ...) {
   #----------------------------------
   # handle arguments
   theModel <- .cvamModel( obj, saturated )
   method <- match.arg(method)
   control <- do.call( "cvamControl", control ) 
   omitData <- as.logical( omitData )[1L]
   stopifnot( inherits(prior, "cvamPrior") )
   #----------------------------------
   # get the input data model frame and frequencies
   mc <- match.call( expand.dots=FALSE )
   mc[[1]] <- quote(stats::model.frame)
   m <- match( c("obj", "data", "freq"), names(mc), nomatch=0L )
   mc <- mc[ c(1L,m) ]
   names(mc)[2] <- "formula"
   mc[[2]] <- theModel 
   mc$na.action <- as.name("na.pass")
   mc$drop.unused.levels <- FALSE
   mf <- eval( mc, parent.frame() )
   if( is.null( mf$`(freq)` ) ) {
      freq <- rep(1L, NROW(mf))
      freqSupplied <- FALSE
   } else {
      freq <- mf$`(freq)`
      mf$`(freq)` <- NULL 
      freqSupplied <- TRUE
   }
   freqInt <- as.integer(freq)
   if( any( freq != freqInt ) ) warning( gettext(
      "Some frequencies changed when integerized" ), domain=NA )
   badV <- names(mf)[ ! unlist( lapply( mf, is.factor ) ) ]
   if( length(badV) > 0L ) stop( gettextf(
      "Variable '%s' is not a factor", badV[1L] ), domain=NA )
   if( any(is.na(freq)) ) stop( gettext(
      "Missing values in 'freq' not allowed"), domain = NA )
   #----------------------------------
   # coerce model frame variables to coarsened and aggregate 
   # to create mfSeen
   mf <- lapply( mf, FUN=coarsened, warnIfCoarsened=FALSE)
#   mf <- lapply( mf, FUN=sticky::sticky )
   if( any( attr(theModel, "fixed") ) ) {
      mfFixed <- mf[ attr(theModel, "fixed") ]
      mfFixed <- lapply( mfFixed, dropCoarseLevels )
      mfFixed <- lapply( mfFixed, is.na )
      mfFixed <- unlist( lapply( mfFixed, any ) )
      if( any(mfFixed) ) { 
         fV <- attr(theModel, "vars")[ attr(theModel, "fixed") ]
         fVbad <- fV[mfFixed]
         stop( gettextf( 
            "Non-modeled variable '%s' contains missing or coarsened values",
            fVbad[1L] ), domain = NA ) 
      }
   }
   mf$freq <- freqInt
   form <-  as.formula( paste( "freq ~", 
      paste( attr(theModel, "vars"), collapse="+" ) ), 
      env = environment(theModel) )
   mfSeen <- aggregate(form, FUN=sum, data=mf)
   if( method == "mfSeen" ) return(mfSeen)
   #----------------------------------
   # information describing the coarsened factors in mfSeen
   freqSeen <- mfSeen$freq
   mfSeen$freq <- NULL
   nBaseLevels <- unlist( lapply( mfSeen, FUN=nBaseLevels) )
   nCoarseLevels <- unlist( lapply( mfSeen, FUN=nCoarseLevels) )
   nLevels <- unlist( lapply( mfSeen, FUN=nlevels) )
   fixed <- attr(theModel, "fixed") 
   stopifnot( all( names(fixed) == names(nLevels) ) )
   latent <- unlist( lapply( mfSeen, FUN=is.latentFactor) )
   #----------------------------------
   # create mfTrue
   mfNoData <- mfSeen[integer(),,drop=FALSE]
   mf0 <- data.frame( lapply( mfNoData, FUN=dropCoarseLevels ) )
   mfContr <- lapply( mf0, FUN=contrasts ) # store the contrasts
   mfTrue <- as.data.frame( xtabs( ~., data=mf0), responseName="freq" )
   mfTrue$freq <- NULL
   for( i in seq_along(mfContr) ){
      vN <- names(mfContr)[i]
      contrasts( mfTrue[[vN]] ) <- mfContr[[vN]]
   } 
   if( method == "mfTrue" ) {
      mfTrue$freq <- as.numeric(NA)
      return(mfTrue)
   }
   #--------------------------------------------
   # create the model matrix and offset, which have zero size if
   # saturated=TRUE
   modelMatrixSupplied <- ! is.null(modelMatrix)
   offsetSupplied <- ! is.null(offset)
   if( attr(theModel, "saturated") ) {
      if( modelMatrixSupplied ) warning( gettextf(
         "modelMatrix ignored because 'saturated=TRUE'" ), domain = NA )
      modelMatrix <- model.matrix( theModel, data=mf0 )  
      if( NCOL(modelMatrix) < NROW(mfTrue) ) warning( gettext(
         "In formula 'obj': model does not appear to be saturated" ),
         domain = NA )
      modelMatrix <- matrix( numeric(), 0L, 0L )
      if( offsetSupplied ) warning( gettextf(
         "offset ignored because 'saturated=TRUE'" ), domain = NA )
      offset <- numeric(0L)   
   } else {
      if( ! modelMatrixSupplied ) {
         modelMatrix <- model.matrix( theModel, data=mf0 )
      } else {
         if( NROW( modelMatrix ) != NROW( mfTrue ) ) stop( gettextf(
            "modelMatrix has incorrect number of rows, should be %i",
            NROW(mfTrue) ), domain = NA )
      }
      if( NCOL(modelMatrix) > control$ncolMaxMM ) stop( gettextf(
         "modelMatrix has %i columns, exceeding 'control$ncolMaxMM'",
         NCOL(modelMatrix) ), domain = NA )
      if( ! modelMatrixSupplied ) modelMatrix <-
         model.matrix( theModel, data=mfTrue )
      if( ! offsetSupplied ) offset <- rep(0, NROW(mfTrue) )
      if( length(offset) != NROW(mfTrue) ) stop( gettextf(
         "offset has incorrect length, should be %i",
         NROW(mfTrue) ), domain = NA ) 
   }
   if( method == "modelMatrix" ) {
      if( attr(theModel, "saturated" ) ) stop( gettext(
         "Model is saturated; there is no modelMatrix" ), domain = NA )
      return( modelMatrix )
   }
   #--------------------------------------------
   # create strZero
   if( is.null(strZero) ) {
      strZero <- rep( FALSE, NROW(mfTrue) )
      strZeroSupplied <- FALSE
   } else {
      strZeroSupplied <- TRUE
   }
   stopifnot( is.logical(strZero) )
   if( length(strZero) != NROW(mfTrue) ) stop( gettextf(
      "strZero has incorrect length, should be %i",
      NROW(mfTrue) ), domain = NA ) 
   if( all( strZero ) ) stop( gettext(
      "All cells are structural zeros, which means nothing is possible"),
      domain = NA )
   #--------------------------------------------
   # ensure modelMatrix has full rank
   if( ! attr(theModel, "saturated") ) {
      if( qr(modelMatrix)$rank < NCOL(modelMatrix) ) stop( gettext(
         "modelMatrix does not have full rank" ), domain=NA )
      if( qr(modelMatrix[(!strZero),,drop=FALSE])$rank <
         NCOL(modelMatrix) ) stop( gettext(
         "modelMatrix does not have full rank due to structural zeros" ), 
         domain=NA )
   }      
   #--------------------------------------------
   # ensure the model contains a constant
   if( ! attr(theModel, "saturated") ) {
      suppressWarnings( fit <- lm(
         rep(1,NROW(modelMatrix)) ~ -1 + modelMatrix,
         subset = ! strZero ) )
      if( sum( fit$residuals^2 ) >= control$critModelCheck ) 
         stop( gettext(
         "In formula 'obj': modelMatrix does not include a constant" ),
         domain = NA )
   }
   #--------------------------------------------
   # check whether model is saturated wrt fixed variables
   if( ! attr(theModel, "saturated") & any( attr(theModel, "fixed") ) ) {
      aV <- attr(theModel, "vars")
      fV <- aV[ attr(theModel, "fixed") ]
      form <-  as.formula( paste( "~", paste( fV, collapse="*" ) ), 
         env = environment(theModel) )
      fM <- model.matrix( form, data=mfTrue )
      suppressWarnings( fit <- 
         lm( fM ~ - 1 + modelMatrix, subset = ! strZero ) )
      if( sum( fit$residuals^2 ) >= control$critModelCheck ) 
         stop( gettext(
         "In formula 'obj': model not saturated wrt fixed variables" ),
         domain = NA )
   }
   #--------------------------------------------
   # process the prior
   pF <- mfSeen[ seq_along(prior$nuggets),,drop=FALSE]
   pS <- .cvamPriorSummaries( prior, pF )
   pF <- pS$pF
   flatten <- pS$flatten
   ridge <- pS$ridge
   priorDataFreq <- pF$freq
   pF$freq <- NULL
   priorDataFreqInt <- pF$freqInt
   pF$freqInt <- NULL
   if( any( priorDataFreq != priorDataFreqInt ) ) message( gettext(
      "Some prior nuggets changed when integerized" ), domain = NA )
   mfPrior <- pF
   mfPrior$freq <- priorDataFreqInt
   if( method == "mfPrior" ) return(mfPrior)
   #--------------------------------------------
   # prob, beta, vhatBeta
   if( is.null(startVal) ) {
      startValSupplied <- FALSE
      probMean <- prob <- numeric( NROW(mfTrue) )
      if( attr(theModel, "saturated") ) {
         betaMean <- beta <- numeric(0L)
         betaCovMat <- vhatBeta <- matrix( numeric(), 0L, 0L )
      } else {
      betaMean <- beta <- structure( numeric( NCOL(modelMatrix) ),
         names = colnames(modelMatrix) )
      betaCovMat <- vhatBeta <- structure(
         matrix( numeric(1L), NCOL(modelMatrix), NCOL(modelMatrix) ),
         rownames = names(beta), colnames = names(beta) )
      }
   } else {
      startValSupplied <- TRUE
      if( attr(theModel, "saturated") ) {
         if( length(startVal) != NROW(mfTrue) ) stop( gettextf(
            "startVal has incorrect length, should be %i", NROW(mfTrue) ),
            domain = NA )
         prob <- as.double( startVal )
         probMean <- numeric( length(prob) )
         betaMean <- beta <- numeric(0L)
         betaMean[] <- as.double(0)
         betaCovMat <- vhatBeta <- matrix( numeric(), 0L, 0L )
      } else {
         if( length(startVal) != NCOL(modelMatrix) ) stop( gettextf(
            "startVal has incorrect length, should be %i", 
            NCOL(modelMatrix) ), domain = NA )
         beta <- as.double( startVal )
         betaMean <- numeric( length(beta) )
         betaCovMat <- vhatBeta <- structure(
            matrix( numeric(1L), NCOL(modelMatrix), NCOL(modelMatrix) ),
            rownames = names(beta), colnames = names(beta) )
         probMean <- prob <- numeric( NROW(mfTrue) )
      }
   }
   #--------------------------------------------
   # process estimates
   argEst <- estimate
   if( is.null(estimate) ) {
      estimate <- list()
   } else {
      if( is.list(estimate)  ){
         ios <- unlist( lapply( estimate, .isOneSided ) )
         if( any( !ios ) ) stop( gettext(
            "Some elements of 'estimate' are not one-sided formulas"),
            domain = NA )
      } else {
         if( ! .isOneSided(estimate) ) stop( gettext(
            "'estimate' is not a one-sided formula"), domain = NA )
         estimate <- list(estimate)
      }
   }
   est <- lapply( estimate, .cvamEstimateRequest, theModel )
   est <- lapply( est, .cvamEstimate, theModel, mfNoData, mf0 )
   estSumm <- .cvamEstimateSummaries(est)
   #--------------------------------------------
   if( method == "MCMC" && control$typeMCMC == "RWM" ) stop( gettext(
      "Cannot run random-walk Metropolis without running EM first"),
      domain = NA )
   if( method == "approxBayes" ) stop( gettext(
      "Cannot run approximate Bayes without running EM first"),
      domain = NA )
   betaHat <- beta
   vhatBetaRWM <- vhatBeta
   #--------------------------------------------
   # prepare objects for .Fortran call
   modelTypeInt <- if( attr(theModel, "saturated" ) ) 1L else 2L
   methodInt <- match( method,  c("EM", "MCMC", "approxBayes" ), nomatch=0L )
   stopifnot( methodInt != 0L )
   packedMap <- unlist( lapply( mfSeen, FUN=.packedMapping ) )
   storage.mode( packedMap ) <- "integer"
   # Note: what is passed to Fortran is actually the seen data, not
   # the data supplied to this function, but for legacy reasons the
   # Fortran code calls it input data
   inputData <- data.matrix(mfSeen)
   storage.mode( inputData ) <- "integer"
   inputDataFreqInt <- mfSeen$freq <- freqSeen
   storage.mode( inputDataFreqInt ) <- "integer" 
   nLevelsMatrix <- cbind( nBaseLevels=nBaseLevels, 
      nCoarseLevels=nCoarseLevels, nLevels=nLevels,
      fixed=as.integer(fixed) ) # 0=FALSE, 1=TRUE
   storage.mode(nLevelsMatrix) <- "integer"
   storage.mode(modelMatrix) <- "double"
   storage.mode(offset) <- "double"
   strZeroInt <- as.integer(strZero)
   storage.mode(strZeroInt) <- "integer"
   storage.mode(flatten) <- "double"
   storage.mode(ridge) <- "double"
   priorData <- data.matrix(pF)
   storage.mode(priorData) <- "integer"
   storage.mode(priorDataFreq) <- "double"
   storage.mode(priorDataFreqInt) <- "integer"
   storage.mode(prob) <- "double"
   storage.mode(beta) <- "double"
   maxIter <- 0L
   if( method == "EM" ) maxIter <- control$iterMaxEM
   if( method == "MCMC" ) maxIter <- control$iterMCMC
   if( method == "approxBayes" ) maxIter <- as.integer( 
      control$iterApproxBayes )
   dimVec <- c(
       nrowInputData = NROW(inputData),    #  1
       ncolInputData = NCOL(inputData),    #  2
       nPackedMap = length(packedMap),     #  3
       nrowTrueData = NROW(mfTrue),        #  4
       nrowMM = NROW(modelMatrix),         #  5
       ncolMM = NCOL(modelMatrix),         #  6
       nrowPriorData = NROW(priorData),    #  7
       maxIter = maxIter )                 #  8
   storage.mode(dimVec) <- "integer"
   ctrlReal <- c( startValJitter = control$startValJitter,   # 1 
      critEM = control$critEM,                               # 2
      critNR = control$critNR,                               # 3
      critBoundary = control$critBoundary )                  # 4
   storage.mode(ctrlReal) <- "double"
   ctrlInt <- c( iterMaxEM = control$iterMaxEM,              # 1
      iterMaxNR = control$iterMaxNR,                         # 2
      startValUseInt = as.integer( startValSupplied ),       # 3
      startValDefaultInt = match( control$startValDefault,   # 4
         c("center", "uniform"), nomatch=0L ),
      excludeAllNAint = as.integer( control$excludeAllNA ) ) # 5
   storage.mode(ctrlInt) <- "integer"
   omitDataInt <- as.integer( omitData )
   nPackedSEs <- if( attr(theModel, "saturated") ) 0L else
      estSumm$nPackedEstimates
   dimVecEst <- c(
      nEstimates = estSumm$nEstimates,                       # 1
      nVarEstimateTot = estSumm$nVarEstimateTot,             # 2
      nPackedEstimates = estSumm$nPackedEstimates,           # 3
      nPackedSEs = nPackedSEs )                              # 4
   storage.mode( dimVecEst ) <- "integer"
   #--------------------------------------------
   # prepare objects to pass to Fortran related to MCMC and approxBayes
   iterMCMC <- as.integer( ceiling( control$iterMCMC / control$thinMCMC ) *
      control$thinMCMC )
   if( method == "MCMC" ) {
      seriesLength <- iterMCMC / control$thinMCMC
      nImpute <- if( control$imputeEvery == 0L ) 0L else
         floor( iterMCMC / control$imputeEvery )
   } else if( method == "approxBayes" ) { 
      seriesLength <- control$iterApproxBayes
      nImpute <- if( control$imputeApproxBayes )
         control$iterApproxBayes else 0L
   } else {
      seriesLength <- 0L
      nImpute <- 0L
   }
   ctrlMCMCInt <- c(
      iterMCMC = iterMCMC,                    # 1
      burnMCMC = control$burnMCMC,            # 2
      thinMCMC = control$thinMCMC,            # 3
      imputeEvery = control$imputeEvery,      # 4
      saveProbSeriesInt = as.integer( control$saveProbSeries ),   # 5
      typeMCMCInt = match( control$typeMCMC,  # 6
         c("DA","RWM"), nomatch=0L ),
      stuckLimit = control$stuckLimit,        # 7
      iterApproxBayes = control$iterApproxBayes, # 8
      imputeApproxBayesInt = as.integer( control$imputeApproxBayes ) ) #9
   storage.mode(ctrlMCMCInt) <- "integer"
   ctrlMCMCReal <- c(
      dfDA = control$tuneDA[1L],             # 1
      stepSizeDA = control$tuneDA[2L],       # 2
      scaleFacDA = control$tuneDA[3L],       # 3
      dfRWM = control$tuneRWM[1L],           # 4
      scaleFacRWM = control$tuneRWM[2L] )    # 5
   storage.mode(ctrlMCMCReal) <- "double"
   if( method == "MCMC" | method == "approxBayes" ) {
      dimVecMCMC <- c(
         seriesLengthBeta = seriesLength,                            # 1
         seriesLengthProb = if( control$saveProbSeries )             # 2
         seriesLength else 0L,
         nImpute = nImpute )                                         # 3
   } else { 
      dimVecMCMC <- c(
         seriesLengthBeta = 0L,                          # 1
         seriesLengthProb = 0L,                          # 2
         nImpute = 0L )                                  # 3
   }
   storage.mode( dimVecMCMC ) <- "integer"
   betaSeries <- matrix( numeric(1L), 
      dimVecMCMC[["seriesLengthBeta"]], length(beta) )
   colnames(betaSeries) <- names(beta)
   probSeries <- matrix( numeric(1L),
      dimVecMCMC[["seriesLengthProb"]], length(prob) )
   imputedFreqInt <- matrix( integer(1L), NROW(mfTrue),
      dimVecMCMC[["nImpute"]] )
   logPSeries <- numeric( dimVecMCMC[["seriesLengthBeta"]] )
   packedEstimatesSeries <- matrix( numeric(1L),
      dimVecMCMC[["seriesLengthBeta"]], estSumm$nPackedEstimates )
   #--------------------------------------------
   # correct problem with length of loglikPadded and logPPadded 
   # related to burn-in
   if( method=="MCMC" ) {
      maxIter <- as.integer(iterMCMC + control$burnMCMC)
      dimVec[["maxIter"]] <- maxIter
   }
   #--------------------------------------------
   # create a matrix for holding message codes
   msg.len.max <- 40L
   msg.codes <- matrix( 0L, msg.len.max, 17L )
   #--------------------------------------------
   tmp <- .Fortran("fit_cvam_model",
      # inputs
      modelTypeInt = modelTypeInt,
      methodInt = methodInt,
      dimVec = dimVec,
      inputData = inputData,
      inputDataFreqInt = inputDataFreqInt,
      nLevelsMatrix = nLevelsMatrix,
      packedMap = packedMap,
      modelMatrix = modelMatrix,
      offset = offset,
      strZeroInt = strZeroInt,
      flatten = flatten, 
      ridge = ridge,
      priorData = priorData,
      priorDataFreqInt = priorDataFreqInt,
      ctrlReal = ctrlReal,
      ctrlInt = ctrlInt,
      omitDataInt = omitDataInt,
      dimVecEst = dimVecEst,
      estimateInfo = estSumm$estimateInfo,
      estimateVarInfo = estSumm$estimateVarInfo,
      ctrlMCMCInt = ctrlMCMCInt,
      ctrlMCMCReal = ctrlMCMCReal,
      dimVecMCMC = dimVecMCMC,
      # inouts
      prob = prob,
      beta = beta,
      betaHat = betaHat,
      vhatBetaRWM = vhatBetaRWM,
      # outputs
      iter = integer(1L),
      convergedInt = integer(1L),
      maxDiff = numeric(1L),
      logliklogPPadded = matrix(numeric(1L), maxIter, 2L),
      lambda = numeric( length(beta) ),
      freq = numeric( NROW(mfTrue) ),
      freqMean = numeric( NROW(mfTrue) ),
      freqInt = integer( NROW(mfTrue) ),
      score = numeric( NCOL(modelMatrix) ),
      vhatBeta = vhatBeta,
      probMean = probMean,
      betaMean = betaMean,
      betaCovMat = betaCovMat,
      totalFreqUsePrior = numeric(1L),
      totalFreqUseDataInt = integer(1L),
      degreesOfFreedom = integer(1L),
      packedEstimates = numeric(estSumm$nPackedEstimates),
      packedEstimatesMean = numeric(estSumm$nPackedEstimates),
      packedSEs = numeric(nPackedSEs),
      betaSeries = betaSeries,
      probSeries = probSeries,
      logPSeries = logPSeries,
      imputedFreqInt = imputedFreqInt,
      packedEstimatesSeries = packedEstimatesSeries,
      nActual = integer(3L),
      mhAcceptRate = numeric(1L),
      startLogP = numeric(1L),
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
   }
   else{
      msg <- paste0( msg.lines, collapse="\n" )
   }
   msg <- paste( msg, "\n", sep="")
   if( msg!= "OK\n" ) cat( paste("Note: ", msg, sep="") )
   #--------------------------------------------
   if( tmp$status != 0 ) stop( gettext( 
      "Procedure aborted" ), domain = NA )
   #--------------------------------------------
    if( ( method == "EM" ) & ( tmp$iter > 0L ) & 
       ( ! tmp$converged ) ) warning( gettextf(
       "Procedure failed to converge by iteration %i", tmp$iter ), 
       domain = NA )
   #--------------------------------------------
   tmp$nIterActual <- tmp$nActual[1L]
   tmp$nSampleActual <- tmp$nActual[2L]
   tmp$nImpActual <- tmp$nActual[3L]
   #--------------------------------------------
   if( length(est) == 0L ) {
      est <- NULL }
   else {
      if( method == "EM" ) {
         for( i in seq_along(est) ) {
            st <- estSumm$estimateInfo[i,3L]
            fin <- estSumm$estimateInfo[i,4L]
            est[[i]]$prob <- tmp$packedEstimates[st:fin]
            if( length(tmp$packedSEs) > 0L ) est[[i]]$SE <- 
               tmp$packedSEs[st:fin]
         }
      } else if( method == "MCMC" || method == "approxBayes" ) {
         if( tmp$nSampleActual <= 0L ) {
	    warning( gettext(
            "No estimates available; insufficient iterations"), domain = NA )
            est <- NULL
         } else {
            for( i in seq_along(est) ) {
               st <- estSumm$estimateInfo[i,3L]
               fin <- estSumm$estimateInfo[i,4L]
               est[[i]]$prob <- apply(  
                  tmp$packedEstimatesSeries[1:tmp$nSampleActual,
                  st:fin, drop=FALSE], 2, mean )
               est[[i]]$SE <- sqrt( apply(  
                  tmp$packedEstimatesSeries[1:tmp$nSampleActual,
                  st:fin, drop=FALSE], 2, var ) )
               est[[i]]$prob.lower <- apply(  
                  tmp$packedEstimatesSeries[1:tmp$nSampleActual,
                  st:fin, drop=FALSE], 2, quantile, (1-control$confidence)/2 )
               est[[i]]$prob.upper <- apply(  
                  tmp$packedEstimatesSeries[1:tmp$nSampleActual,
                  st:fin, drop=FALSE], 2, quantile, 1-(control$confidence/2) )
            }
         }
      }
   }
   if( ! is.null(est) ) {
      nDraws <- if( method == "MCMC" | method == "approxBayes") 
         tmp$nSampleActual else 0L
      est <- lapply( est, .formatEstimate, control$confidence,
         control$probRound, control$probDigits,
         fromSeries = ( method == "MCMC" | method == "approxBayes" ),
         fromMCMC = ( method == "MCMC" | method == "approxBayes" ),
            nDraws = nDraws, meanSeries=TRUE )
      est <- if( length(est) == 1L ) est[[1L]] else est
      est <- .cvamEstimateList( est )
   }
   #--------------------------------------------
   tmp$beganAtMode <- FALSE
   tmp$atMode <- ( method == "EM" ) & as.logical( tmp$convergedInt )
   tmp$converged <- if( method == "EM")
      as.logical( tmp$convergedInt ) else NULL
   loglikPadded <- tmp$logliklogPPadded[,1L]
   logPPadded <- tmp$logliklogPPadded[,2L]
   tmp$loglik <- if( tmp$iter > 0L ) 
      loglikPadded[1:tmp$iter] else numeric()
   tmp$logP <- if( tmp$iter > 0L ) 
      logPPadded[1:tmp$iter] else numeric()
   mfTrue$freq <- if( method == "EM" ) tmp$freq else tmp$freqMean
   tmp$mfTrue <- mfTrue
   tmp$mfSeen <- mfSeen
   tmp$mfPrior <- mfPrior
   tmp$model <- theModel
   tmp$method <- method
   tmp$prior <- prior
   tmp$control <- control
   tmp$omitData <- omitData
   tmp$startValSupplied <- startValSupplied
   tmp$offsetSupplied <- offsetSupplied
   tmp$freqSupplied <- freqSupplied
   tmp$strZeroSupplied <- strZeroSupplied
   tmp$modelMatrixSupplied <- modelMatrixSupplied
   tmp$strZero <- strZero
   tmp$argEst <- argEst
   tmp$estimates <- est
   tmp$latent <- latent
   tmp$call <- match.call()
   #--------------------------------------------
   tmp$betaHat <- if( tmp$atMode ) tmp$beta else NULL
   tmp$vhatBetaRWM <- if( tmp$atMode ) tmp$vhatBeta else NULL
   #--------------------------------------------
   if( ( method == "MCMC" || method == "approxBayes" ) && 
      ( tmp$nSampleActual >= 1L ) ){
      tmp$logPseries <- tmp$logPseries[ 1:tmp$nSampleActual ]
   } else {
      tmp$logPseries <- NULL
   }
   if( ( method == "MCMC" || method == "approxBayes" ) &&
      ( tmp$nSampleActual >= 1L ) && 
      ( ! attr(theModel, "saturated") ) ){
      tmp$betaSeries <- tmp$betaSeries[ 1:tmp$nSampleActual,,drop=FALSE ]
   } else {
      tmp$betaSeries <- NULL
   }
   if( ( method == "MCMC" || method == "approxBayes" ) &&
      ( tmp$nSampleActual >= 1L ) &&
      ( control$saveProbSeries ) ){
      tmp$probSeries <- tmp$probSeries[ 1:tmp$nSampleActual,,drop=FALSE ]
   } else {
      tmp$probSeries <- NULL
   }
   if( ( method == "MCMC" || method == "approxBayes" ) && 
      ( tmp$nSampleActual <= 0L ) ){
      tmp$betaMean <- tmp$betaCovMat <- NULL
   }
   if( ( method == "MCMC" || method == "approxBayes" ) && 
      ( tmp$nSampleActual >= 1L ) &&
      ( !is.null(est) ) ){
      tmp$packedEstimatesSeries <-
         tmp$packedEstimatesSeries[1:tmp$nSampleActual,,drop=FALSE]
   } else {
      tmp$packedEstimatesSeries <- NULL
   }
   #--------------------------------------------
   structure( tmp, class = c("cvam", "list") )
}



cvam.cvam <- function(obj, method=obj$method, control=NULL, startVal=NULL,
   estimate=NULL, ...) {
   #----------------------------------
   stopifnot( inherits(obj, "cvam") )
   # disable the possibility to change
   # prior and omitData
   prior <- obj$prior
   omitData <- obj$omitData
   stopifnot( method %in% c("EM", "MCMC", "approxBayes", "mfSeen", "mfTrue",
      "mfPrior", "modelMatrix") )
   omitData <- as.logical(omitData)[1L]
   stopifnot( inherits(prior, "cvamPrior") )
   #----------------------------------
   if( method == "modelMatrix" ) {
      if( attr(obj$model, "saturated" ) ) stop( gettext(
         "Model is saturated; there is no modelMatrix" ), domain = NA )
      return( obj$modelMatrix )
   }
   #----------------------------------
   if( is.null(control) ) {
      control <- obj$control
   } else {
      stopifnot( is.list(control) )
      ctrl <- obj$control
      for( i in seq_along(control) ) ctrl[[ names(control)[i] ]] <- 
         control[[i]]
      control <- do.call( "cvamControl", ctrl )
   }
   #----------------------------------
   prob <- obj$prob
   beta <- obj$beta
   if( ! is.null(startVal) ) {
      if( attr(obj$model, "saturated") ) {
         if( length(startVal) != NROW(obj$mfTrue) ) stop( gettextf(
            "startVal has incorrect length, should be %i", NROW(obj$mfTrue) ),
            domain = NA )
         prob[] <- as.double( startVal )
      } else {
         if( length(startVal) != NCOL(obj$modelMatrix) ) stop( gettextf(
            "startVal has incorrect length, should be %i", 
            NCOL(obj$modelMatrix) ), domain = NA )
         beta[] <- as.double( startVal )
      }
   }
   doingRWM <-  ( ! attr(obj$model, "saturated") ) &
      method == "MCMC" & control$typeMCMC == "RWM" 
   if( doingRWM & is.null(obj$vhatBetaRWM) )
      stop( gettext(
      "Cannot run random-walk Metropolis without first getting a mode from EM"),
          domain = NA )
   doingApproxBayes <- ( ! attr(obj$model, "saturated") ) &
      method == "approxBayes"
   if( doingApproxBayes & ( is.null(obj$vhatBetaRWM) | is.null(obj$betaHat) ) )
      stop( gettext(
      "Cannot run approximate Bayes without first getting a mode from EM"),
          domain = NA )
   if( attr(obj$model, "saturated") ) {
      betaHat <- numeric(0L)
      vhatBeta <- vhatBetaRWM <- matrix( numeric(), 0L, 0L )
   } else {
      if( doingRWM | doingApproxBayes ) {
         betaHat <- obj$betaHat
         vhatBeta <- vhatBetaRWM <- obj$vhatBetaRWM
         vhatBeta[] <- as.double(0)
      } else {
         betaHat <- numeric( NCOL(obj$modelMatrix) )
         vhatBeta <- vhatBetaRWM <- structure(
            matrix( numeric(1L), NCOL(obj$modelMatrix), NCOL(obj$modelMatrix) ),
            rownames = names(beta), colnames = names(beta) )
      }
   }
   #
   startValSupplied <- TRUE
   probMean <- prob
   probMean[] <- as.double(0)
   betaMean <- beta
   betaMean[] <- as.double(0)
   betaCovMat <- vhatBeta
   #----------------------------------
   mfSeen <- obj$mfSeen
   if( method == "mfSeen" ) return(mfSeen)
   freqSeen <- mfSeen$freq
   mfSeen$freq <- NULL
   mfNoData <- mfSeen[integer(),,drop=FALSE]
   mf0 <- data.frame( lapply( mfNoData, FUN=dropCoarseLevels ) )
   mfContr <- lapply( mf0, FUN=contrasts ) # store the contrasts
   mfTrue <- as.data.frame( xtabs( ~., data=mf0), responseName="freq" )
   mfTrue$freq <- NULL
   for( i in seq_along(mfContr) ){
      vN <- names(mfContr)[i]
      contrasts( mfTrue[[vN]] ) <- mfContr[[vN]]
   } 
   if( method == "mfTrue" ) {
      mfTrue$freq <- as.numeric(NA)
      return(mfTrue)
   }
   #----------------------------------
   pF <- mfSeen[ seq_along(prior$nuggets),,drop=FALSE]
   pS <- .cvamPriorSummaries( prior, pF )
   pF <- pS$pF
   flatten <- pS$flatten
   ridge <- pS$ridge
   priorDataFreq <- pF$freq
   pF$freq <- NULL
   priorDataFreqInt <- pF$freqInt
   pF$freqInt <- NULL
   if( any( priorDataFreq != priorDataFreqInt ) ) message( gettext(
      "Some prior nuggets changed when integerized" ), domain = NA )
   mfPrior <- pF
   mfPrior$freq <- priorDataFreqInt
   if( method == "mfPrior" ) return(mfPrior)
   #----------------------------------
   argEst <- estimate
   if( is.null(estimate) ) estimate <- argEst <- obj$argEst
   if( is.null(estimate) ) {
      estimate <- list()
   } else {
      if( length(estimate) > 0L ) {
         if( is.list(estimate) ){
            ios <- unlist( lapply( estimate, .isOneSided ) )
            if( any( !ios ) ) stop( gettext(
               "Some elements of 'estimate' are not one-sided formulas"),
               domain = NA )
         } else {
            if( ! .isOneSided(estimate) ) stop( gettext(
               "'estimate' is not a one-sided formula"), domain = NA )
            estimate <- list(estimate)
         }
      }
   }
   est <- lapply( estimate, .cvamEstimateRequest, obj$model )
   est <- lapply( est, .cvamEstimate, obj$model, mfNoData, mf0 )
   estSumm <- .cvamEstimateSummaries(est)
   #--------------------------------------------
   # prepare objects for .Fortran call
   methodInt <- match( method,  c("EM", "MCMC", "approxBayes" ), nomatch=0L )
   stopifnot( methodInt != 0L )
   storage.mode(flatten) <- "double"
   storage.mode(ridge) <- "double"
   priorData <- data.matrix(pF)
   storage.mode(priorData) <- "integer"
   storage.mode(priorDataFreq) <- "double"
   storage.mode(priorDataFreqInt) <- "integer"
   storage.mode(prob) <- "double"
   storage.mode(beta) <- "double"
   maxIter <- 0L
   if( method == "EM" ) maxIter <- control$iterMaxEM
   if( method == "MCMC" ) maxIter <- control$iterMCMC
   if( method == "approxBayes" ) maxIter <- as.integer( 
      control$iterApproxBayes )
   dimVec <- c(
       nrowInputData = NROW(obj$inputData),    #  1
       ncolInputData = NCOL(obj$inputData),    #  2
       nPackedMap = length(obj$packedMap),     #  3
       nrowTrueData = NROW(mfTrue),            #  4
       nrowMM = NROW(obj$modelMatrix),         #  5
       ncolMM = NCOL(obj$modelMatrix),         #  6
       nrowPriorData = NROW(priorData),        #  7
       maxIter = maxIter )                     #  8
   storage.mode(dimVec) <- "integer"
   ctrlReal <- c( startValJitter = control$startValJitter,   # 1 
      critEM = control$critEM,                               # 2
      critNR = control$critNR,                               # 3
      critBoundary = control$critBoundary )                  # 4
   storage.mode(ctrlReal) <- "double"
   ctrlInt <- c( iterMaxEM = control$iterMaxEM,              # 1
      iterMaxNR = control$iterMaxNR,                         # 2
      startValUseInt = as.integer( startValSupplied ),       # 3
      startValDefaultInt = match( control$startValDefault,   # 4
         c("center", "uniform"), nomatch=0L ),
      excludeAllNAint = as.integer( control$excludeAllNA ) ) # 5
   storage.mode(ctrlInt) <- "integer"
   omitDataInt <- as.integer( omitData )
   nPackedSEs <- if( attr(obj$model, "saturated") ) 0L else
      estSumm$nPackedEstimates
   dimVecEst <- c(
      nEstimates = estSumm$nEstimates,                       # 1
      nVarEstimateTot = estSumm$nVarEstimateTot,             # 2
      nPackedEstimates = estSumm$nPackedEstimates,           # 3
      nPackedSEs = nPackedSEs )                              # 4
   storage.mode( dimVecEst ) <- "integer"
   #--------------------------------------------
   # prepare objects to pass to Fortran related to MCMC and approxBayes
   iterMCMC <- as.integer( ceiling( control$iterMCMC / control$thinMCMC ) *
      control$thinMCMC )
   if( method == "MCMC" ) {
      seriesLength <- iterMCMC / control$thinMCMC
      nImpute <- if( control$imputeEvery == 0L ) 0L else
         floor( iterMCMC / control$imputeEvery )
   } else if( method == "approxBayes" ) { 
      seriesLength <- control$iterApproxBayes
      nImpute <- if( control$imputeApproxBayes )
         control$iterApproxBayes else 0L
   } else {
      seriesLength <- 0L
      nImpute <- 0L
   }
   ctrlMCMCInt <- c(
      iterMCMC = iterMCMC,                    # 1
      burnMCMC = control$burnMCMC,            # 2
      thinMCMC = control$thinMCMC,            # 3
      imputeEvery = control$imputeEvery,      # 4
      saveProbSeriesInt = as.integer( control$saveProbSeries ),   # 5
      typeMCMCInt = match( control$typeMCMC,  # 6
         c("DA","RWM"), nomatch=0L ),
      stuckLimit = control$stuckLimit,        # 7
      iterApproxBayes = control$iterApproxBayes, # 8
      imputeApproxBayesInt = as.integer( control$imputeApproxBayes ) ) #9
   storage.mode(ctrlMCMCInt) <- "integer"
   ctrlMCMCReal <- c(
      dfDA = control$tuneDA[1L],             # 1
      stepSizeDA = control$tuneDA[2L],       # 2
      scaleFacDA = control$tuneDA[3L],       # 3
      dfRWM = control$tuneRWM[1L],           # 4
      scaleFacRWM = control$tuneRWM[2L] )    # 5
   storage.mode(ctrlMCMCReal) <- "double"
   if( method == "MCMC" | method == "approxBayes" ) {
      dimVecMCMC <- c(
         seriesLengthBeta = seriesLength,                            # 1
         seriesLengthProb = if( control$saveProbSeries )             # 2
         seriesLength else 0L,
         nImpute = nImpute )                                         # 3
   } else { 
      dimVecMCMC <- c(
         seriesLengthBeta = 0L,                          # 1
         seriesLengthProb = 0L,                          # 2
         nImpute = 0L )                                  # 3
   }
   storage.mode( dimVecMCMC ) <- "integer"
   betaSeries <- matrix( numeric(1L), 
      dimVecMCMC[["seriesLengthBeta"]], length(beta) )
   colnames(betaSeries) <- names(beta)
   probSeries <- matrix( numeric(1L),
      dimVecMCMC[["seriesLengthProb"]], length(prob) )
   imputedFreqInt <- matrix( integer(1L), NROW(mfTrue),
      dimVecMCMC[["nImpute"]] )
   logPSeries <- numeric( dimVecMCMC[["seriesLengthBeta"]] )
   packedEstimatesSeries <- matrix( numeric(1L),
      dimVecMCMC[["seriesLengthBeta"]], estSumm$nPackedEstimates )
   #--------------------------------------------
   # correct problem with length of loglikPadded and logPPadded 
   # related to burn-in
   if( method=="MCMC" ) {
      maxIter <- as.integer(iterMCMC + control$burnMCMC)
      dimVec[["maxIter"]] <- maxIter
   }
   #--------------------------------------------
   # create a matrix for holding message codes
   msg.len.max <- 40L
   msg.codes <- matrix( 0L, msg.len.max, 17L )
   #--------------------------------------------
   tmp <- .Fortran("fit_cvam_model",
      # inputs
      modelTypeInt = obj$modelTypeInt,
      methodInt = methodInt,
      dimVec = dimVec,
      inputData = obj$inputData,
      inputDataFreqInt = obj$inputDataFreqInt,
      nLevelsMatrix = obj$nLevelsMatrix,
      packedMap = obj$packedMap,
      modelMatrix = obj$modelMatrix,
      offset = obj$offset,
      strZeroInt = obj$strZeroInt,
      flatten = flatten, 
      ridge = ridge,
      priorData = priorData,
      priorDataFreqInt = priorDataFreqInt,
      ctrlReal = ctrlReal,
      ctrlInt = ctrlInt,
      omitDataInt = omitDataInt,
      dimVecEst = dimVecEst,
      estimateInfo = estSumm$estimateInfo,
      estimateVarInfo = estSumm$estimateVarInfo,
      ctrlMCMCInt = ctrlMCMCInt,
      ctrlMCMCReal = ctrlMCMCReal,
      dimVecMCMC = dimVecMCMC,
      # inouts
      prob = prob,
      beta = beta,
      betaHat = betaHat,
      vhatBetaRWM = vhatBetaRWM,
      # outputs
      iter = integer(1L),
      convergedInt = integer(1L),
      maxDiff = numeric(1L),
      logliklogPPadded = matrix(numeric(1L), maxIter, 2L),
      lambda = numeric( length(beta) ),
      freq = numeric( NROW(mfTrue) ),
      freqMean = numeric( NROW(mfTrue) ),
      freqInt = integer( NROW(mfTrue) ),
      score = numeric( NCOL(obj$modelMatrix) ),
      vhatBeta = vhatBeta,
      probMean = probMean,
      betaMean = betaMean,
      betaCovMat = betaCovMat,
      totalFreqUsePrior = numeric(1L),
      totalFreqUseDataInt = integer(1L),
      degreesOfFreedom = integer(1L),
      packedEstimates = numeric(estSumm$nPackedEstimates),
      packedEstimatesMean = numeric(estSumm$nPackedEstimates),
      packedSEs = numeric(nPackedSEs),
      betaSeries = betaSeries,
      probSeries = probSeries,
      logPSeries = logPSeries,
      imputedFreqInt = imputedFreqInt,
      packedEstimatesSeries = packedEstimatesSeries,
      nActual = integer(3L),
      mhAcceptRate = numeric(1L),
      startLogP = numeric(1L),
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
   }
   else{
      msg <- paste0( msg.lines, collapse="\n" )
   }
   msg <- paste( msg, "\n", sep="")
   if( msg!= "OK\n" ) cat( paste("Note: ", msg, sep="") )
   #--------------------------------------------
   if( tmp$status != 0 ) stop( gettext( 
      "Procedure aborted" ), domain = NA )
   #--------------------------------------------
    if( ( method == "EM" ) & ( tmp$iter > 0L ) & 
       ( ! tmp$converged ) ) warning( gettextf(
       "Procedure failed to converge by iteration %i", tmp$iter ), 
       domain = NA )
   #--------------------------------------------
   tmp$nIterActual <- tmp$nActual[1L]
   tmp$nSampleActual <- tmp$nActual[2L]
   tmp$nImpActual <- tmp$nActual[3L]
   #--------------------------------------------
   if( length(est) == 0L ) {
      est <- NULL }
   else {
      if( method == "EM" ) {
         for( i in seq_along(est) ) {
            st <- estSumm$estimateInfo[i,3L]
            fin <- estSumm$estimateInfo[i,4L]
            est[[i]]$prob <- tmp$packedEstimates[st:fin]
            if( length(tmp$packedSEs) > 0L ) est[[i]]$SE <- 
               tmp$packedSEs[st:fin]
         }
      } else if( method == "MCMC" || method == "approxBayes" ) {
         if( tmp$nSampleActual < 1L ) {
	    warning( gettext(
            "No estimates available; insufficient iterations"), domain = NA )
            est <- NULL
         } else {
            for( i in seq_along(est) ) {
               st <- estSumm$estimateInfo[i,3L]
               fin <- estSumm$estimateInfo[i,4L]
               est[[i]]$prob <- apply(  
                  tmp$packedEstimatesSeries[1:tmp$nSampleActual,
                  st:fin, drop=FALSE], 2, mean )
               est[[i]]$SE <- sqrt( apply(  
                  tmp$packedEstimatesSeries[1:tmp$nSampleActual,
                  st:fin, drop=FALSE], 2, var ) )
               est[[i]]$prob.lower <- apply(  
                  tmp$packedEstimatesSeries[1:tmp$nSampleActual,
                  st:fin, drop=FALSE], 2, quantile, (1-control$confidence)/2 )
               est[[i]]$prob.upper <- apply(  
                  tmp$packedEstimatesSeries[1:tmp$nSampleActual,
                  st:fin, drop=FALSE], 2, quantile, 1-(control$confidence/2) )
            }
         }
      }
   }
   if( ! is.null(est) ) {
      nDraws <- if( method == "MCMC" | method == "approxBayes" ) 
          tmp$nSampleActual else 0L
      est <- lapply( est, .formatEstimate, control$confidence,
         control$probRound, control$probDigits,
         fromSeries = ( method == "MCMC" | method == "approxBayes" ),
         fromMCMC = ( method == "MCMC" | method == "approxBayes" ),
            nDraws = nDraws, meanSeries=TRUE )
      est <- if( length(est) == 1L ) est[[1L]] else est
      est <- .cvamEstimateList( est )
   }
   #--------------------------------------------
   tmp$beganAtMode <- obj$atMode
   tmp$atMode <- ( method == "EM" ) & as.logical( tmp$convergedInt )
   tmp$converged <- if( method == "EM")
      as.logical( tmp$convergedInt ) else NULL
   loglikPadded <- tmp$logliklogPPadded[,1L]
   logPPadded <- tmp$logliklogPPadded[,2L]
   tmp$loglik <- if( tmp$iter > 0L ) 
      loglikPadded[1:tmp$iter] else numeric()
   tmp$logP <- if( tmp$iter > 0L ) 
      logPPadded[1:tmp$iter] else numeric()
   mfTrue$freq <- if( method == "EM" ) tmp$freq else tmp$freqMean
   tmp$mfTrue <- mfTrue
   mfSeen$freq <- freqSeen
   tmp$mfSeen <- mfSeen
   tmp$mfPrior <- mfPrior
   tmp$model <- obj$model
   tmp$method <- method
   tmp$prior <- prior
   tmp$control <- control
   tmp$omitData <- omitData
   tmp$startValSupplied <- startValSupplied
   tmp$offsetSupplied <- obj$offsetSupplied
   tmp$freqSupplied <- obj$freqSupplied
   tmp$strZeroSupplied <- obj$strZeroSupplied
   tmp$modelMatrixSupplied <- obj$modelMatrixSupplied
   tmp$strZero <- obj$strZero
   tmp$estimates <- est
   tmp$argEst <- argEst
   tmp$latent <- obj$latent
   tmp$call <- match.call()
   #--------------------------------------------
   tmp$betaHat <- if( tmp$atMode ) tmp$beta else obj$betaHat
   tmp$vhatBetaRWM <- if( tmp$atMode ) tmp$vhatBeta else obj$vhatBetaRWM
   #--------------------------------------------
   if( ( method == "MCMC" || method == "approxBayes" ) && 
      ( tmp$nSampleActual >= 1L ) ){
      tmp$logPseries <- tmp$logPseries[ 1:tmp$nSampleActual ]
   } else {
      tmp$logPseries <- NULL
   }
   if( ( method == "MCMC" || method == "approxBayes" ) &&
      ( tmp$nSampleActual >= 1L ) &&
      ( ! attr(obj$model, "saturated") ) ){
      tmp$betaSeries <- tmp$betaSeries[ 1:tmp$nSampleActual,,drop=FALSE ]
   } else {
      tmp$betaSeries <- NULL
   }
   if( ( method == "MCMC" || method == "approxBayes" ) &&
      ( tmp$nSampleActual >= 1L ) &&
      ( control$saveProbSeries ) ){
      tmp$probSeries <- tmp$probSeries[ 1:tmp$nSampleActual,,drop=FALSE ]
   } else {
      tmp$probSeries <- NULL
   }
   if( ( method == "MCMC" || method == "approxBayes" ) &&
      ( tmp$nSampleActual <= control$burnMCMC ) ){
      tmp$betaMean <- tmp$betaCovMat <- NULL
   }
   if( ( method == "MCMC" || method == "approxBayes" ) &&
      ( tmp$nSampleActual >= 1L ) &&
      ( !is.null(est) ) ){
      tmp$packedEstimatesSeries <-
         tmp$packedEstimatesSeries[1:tmp$nSampleActual,,drop=FALSE]
   } else {
      tmp$packedEstimatesSeries <- NULL
   }
   #--------------------------------------------
   structure( tmp, class = c("cvam", "list") )
}

print.cvam <- function(x, ...) {
   stopifnot( inherits(x, "cvam") )
   cat( gettext( 'Object of class "cvam", use "summary"' ), sep="\n")
   invisible()
}

summary.cvam <- function(object, showCoef=TRUE, showEstimates=TRUE,
   digits=4L, ...) {
   stopifnot( inherits(object, "cvam") )
   if( inherits(object, "summary.cvam") ) return(object)
   result <- list()
   result$model <- object$model
   result$prior <- object$prior
   result$method <- object$method
   result$control <- object$control
   result$logP <- if( object$iter > 0L )
      object$logP[object$iter] else numeric()
   result$loglik <- if( object$iter > 0L )
      object$loglik[object$iter] else numeric()
   result$totalSampSize <- sum( object$inputDataFreqInt )
   result$priorSampSize <- object$flatten + sum( object$priorDataFreqInt )
   result$nCells <- prod( object$nLevelsMatrix[,1L] )
   result$anyLatent <- any( object$latent )
   result$nCellsNonLatent <-
      prod( object$nLevelsMatrix[!object$latent,1L] )
   result$nZero <- sum( object$strZero )
   result$nDataPatt <- object$dimVec[1L]
   result$degreesOfFreedom <-  if( result$anyLatent &&
     ! attr(object$model, "saturated") )
     object$degreesOfFreedom -
     ( result$nCells - result$nCellsNonLatent ) else
      object$degreesOfFreedom
   result$nParams <- if( attr(object$model, "saturated") )
      result$nCells - result$nZero else NCOL( object$modelMatrix )
   result$totalNSupplied <- sum( object$inputDataFreqInt )
   result$totalN <- object$totalFreqUseDataInt
   result$priorESS <- object$totalFreqUsePrior
   result$startValSupplied <- if( object$method != "approxBayes" )
      object$startValSupplied else NULL
   result$offsetSupplied <- object$offsetSupplied
   if( object$method == "EM" ) {
      result$converged <- object$converged
      result$iter <- object$iter
      result$lenGrad <- if( length(object$score) > 0L ) 
         sqrt( sum(object$score^2) ) else NULL
   } else if( object$method == "MCMC" ) {
      result$iter <- object$nIterActual
      result$discarded <- if( result$iter >= object$control$burnMCMC ) 
         object$control$burnMCMC else ( result$iter - object$control$burnMCMC )
      result$afterBurnIn <- max( result$iter - object$control$burnMCMC, 0 )
      result$thin <- object$control$thinMCMC
      result$imputeEvery <- object$control$imputeEvery
      result$nSampleActual <- object$nSampleActual
      result$nImpActual <- object$nImpActual
      result$probSeriesPresent <- object$control$saveProbSeries
      result$mhAcceptRate <- object$mhAcceptRate
   } else if( object$method == "approxBayes" ) {
      result$iter <- object$nIterActual
      result$nSampleActual <- object$nSampleActual
      result$nImpActual <- object$nImpActual
      result$probSeriesPresent <- object$control$saveProbSeries
   }
   result$coef <- get.coef( object, withSE = TRUE, msgIfNone = FALSE )
   result$showCoef <- if( is.null(result$coef) ) FALSE else showCoef
   #
   result$beganAtMode <- object$beganAtMode
   result$atMode <- object$atMode
   #
   result$impPresent <- FALSE
   if( ( object$method == "MCMC" || object$method == "approxBayes" ) &&
      ( object$nImpActual > 0L ) ) result$impPresent <- TRUE
   result$estimates <- object$estimates
   result$showEstimates <- showEstimates
   result$digits <- digits
   structure( result,
      class = "summary.cvam" )
}

print.summary.cvam <- function(x, ...) {
   stopifnot( inherits(x, "summary.cvam") )
   cat( attr(x$model, "formulaStr"), sep="\n" )
   if( attr(x$model, "saturated") )
      cat( "   fit as a saturated model", sep="\n")
   if( x$offsetSupplied ) cat( "(offset presnt)", sep="\n" )
   cat("\n")
   cat("Prior:", sep="\n")
   print.summary.cvamPrior( summary(x$prior), showNuggets=FALSE, ...)
   if( attr(x$model, "saturated") && ( x$prior$ridge > 0 ) ) {
      cat("Ridge factor was ignored because saturated=TRUE", sep="\n")
   }
   cat("\n")
   cat("Sample size:", sep="\n")
   strA <- format( c( "total N in supplied data =",
      "N from supplied data used in model fit =",
      "prior effective sample size ="), justify="right" )
   strB <- format( c( x$totalNSupplied, x$totalN, x$priorESS),
      justify = "left" )
   cat( paste(strA, strB), sep="\n")
   cat("\n")
   cat("Degrees of freedom:", sep="\n")
   strA <- format( c("patterns of coarsened data =",
      "cells in complete-data table =",
      "cells without latent variables =",
      "structural zero cells =",
      "parameters in Poisson model =",
      "df ="), justify="right" )
   strB <- format( c(x$nDataPatt, x$nCells, x$nCellsNonLatent,
      x$nZero, x$nParams, x$degreesOfFreedom),
      justify = "left" )
   cat( paste(strA, strB), sep="\n")
   cat("\n")
   if( x$method != "approxBayes" ) {
      cat( "Starting values:", sep="\n")
      if( x$startValSupplied ) {
         cat("supplied by user", "\n")
      } else {
         cat( gettextf("default, %s", x$control$startValDefault), sep="\n")
      }
      cat( gettextf("jitter SD = %f", x$control$startValJitter), sep="\n")
      cat("\n")
   }
   if( x$method == "EM" &  x$iter > 0  ) {
      cat( "EM algorithm:", sep="\n")
       if( x$converged ) {
          cat( gettextf( "Converged at iteration %i", x$iter ), sep="\n")
       } else {
          cat( gettextf( "Failed to converge by iteration %i", 
             x$iter ), sep="\n")
       }
       if( !is.null(x$lenGrad) ) 
           cat( gettextf( "Gradient length = %f", x$lenGrad ), sep="\n" )
       cat("\n")
       strA <- format( c("Final logP =", "Final loglik ="), justify="right")
       strB <- format( c(x$logP, x$loglik), justify="left")
       cat( paste(strA, strB), sep="\n")
       cat("\n")
   } else if( x$method == "MCMC" & x$iter > 0 ){
      if( attr( x$model, "saturated" ) ) {
         pStr <- "MCMC: Data augmentation (DA) for saturated model"
         cat( pStr, sep="\n")
      } else {
         if( x$control$typeMCMC == "DA" ) {
            pStr <- "MCMC: Data augumentation (DA) with Metropolis-Hastings"
            cat( pStr, sep="\n")
            cat("\n")
            cat("Tuning parameters:", sep="\n")
            strA <- format( c("proposal df =", 
               "step size =", "scale factor ="), justify="right")
            strB <- c( format( x$control$tuneDA[1L] ), 
               format( x$control$tuneDA[-1L] ) )
            strB <- format( strB, justify="left")
            cat( paste(strA, strB), sep="\n")
	    cat("\n")
         } else if( x$control$typeMCMC == "RWM" ) {
            pStr <- "MCMC: Random-walk Metropolis"
            cat( pStr, sep="\n")
            cat("\n")
            cat("Tuning parameters:", sep="\n")
            strA <- format( c("proposal df =", 
               "scale factor ="), justify="right")
            strB <- c( format( x$control$tuneRWM[1L] ), 
               format( x$control$tuneRWM[2L] ) )
            strB <- format( strB, justify="left")
            cat( paste(strA, strB), sep="\n")
            cat("\n")
         }
      }
      strA <- "Accept rate ="
      strB <- format(x$mhAcceptRate)
      cat( paste(strA, strB), sep="\n")
      cat("\n")
      strA <- format( c("Iterations performed =", 
         "Iterations discarded as burn-in =", 
         "Iterations after burn-in =", 
         "Thinning interval for saved series =",
         "Samples in saved series =",
         "Imputation interval =",
         "Number of imputations stored ="), justify="right")
      strB <- format( c( x$iter,
         x$discarded,
         x$afterBurnIn,
         x$thin,
         x$nSampleActual,
         x$imputeEvery,
         x$nImpActual), justify="left")
      cat( paste(strA, strB), sep="\n")
      cat("\n")
   } else if( x$method == "approxBayes" & x$iter > 0 ){
      pStr <-
         "Approximate Bayes: Independent draws from approximate posterior"
      cat( pStr, sep="\n")
      strA <- format( c("Iterations performed =", 
         "Samples in saved series =",
         "Number of imputations stored ="), justify="right")
      strB <- format( c( x$iter,
         x$nSampleActual,
         x$nImpActual), justify="left")
      cat( paste(strA, strB), sep="\n")
      cat("\n")
   }
   if( x$showCoef ) {
      cat( attr(x$coef, "header"), sep="\n")
      print(x$coef, digits=x$digits, ...)
      cat("\n")
   }
   if( x$showEstimates & ! is.null(x$estimates) ) {
      print(x$estimates, showHeader = TRUE, digits=x$digits, ...)
   }
   invisible()
}

.cvamEstimateRequest <- function( obj, model ) {
   aN <- setdiff( all.names(obj), c("~", "(", "+", "|") )
   aV <- all.vars(obj)
   badN <- setdiff( aN, aV )
   if( length(badN) > 0L ) stop( gettextf(
      "In formula '%s': symbol '%s' not allowed", 
         deparse(obj), badN[1L] ), domain = NA )
   Form <- Formula::Formula(obj)
   nParts <- length( attr(Form, "rhs") )
   if( nParts > 2L ) stop( gettextf( 
     "Multiple '|' found in '%s'", deparse(obj) ), domain = NA )
   randVars <- all.vars( formula(Form, rhs=1L ) )
   fixedVars <- if( nParts == 2L ) 
      all.vars( formula(Form, rhs=2L) )  else character()
   if( length(randVars) == 0L ) stop( gettextf(
      "In formula '%s': nothing to estimate",
      deparse(obj)), domain = NA )
   inBoth <- intersect( randVars, fixedVars )
   if( length(inBoth) > 0L ) stop( gettextf(
      "In formula '%s': variable '%s' appears on both sides of '|'", 
      deparse(obj), inBoth[1L] ), domain = NA )
   allVars <- c( randVars, fixedVars )
   notInModel <- setdiff( allVars, attr(model, "vars") )
   if( length(notInModel) > 0L ) stop( gettextf(
      "In formula '%s': variable '%s' does not appear in the model",
      deparse(obj),  notInModel[1L] ), domain = NA )
   fInModel <- attr(model, "vars")[ attr(model, "fixed") ]
   notModeled <- intersect( randVars, fInModel )
   if( length(notModeled) > 0L ) stop( gettextf(
      "In formula '%s': variable '%s' is fixed in the model",
      deparse(obj),  notModeled[1L] ), domain = NA )
   missingFixedVars <- setdiff( fInModel, fixedVars )
   if( length(missingFixedVars) > 0L ) stop( gettextf(
      "In formula '%s': variable '%s' does not appear after '|'",
      deparse(obj),  missingFixedVars[1L] ), domain = NA )
   formulaStr <- deparse(obj)
   formulaStr <- paste("~ ", substring(formulaStr, 2L), sep="" )
   structure( obj,
      vars = allVars,
      fixed = allVars %in% fixedVars,
      formulaStr = formulaStr,
      .Environment = attr(model, ".Environment"),
      class = c( ".cvamEstimateRequest", class(obj) ) )
}

.cvamEstimate <- function( obj, model, mfNoData, mf0 ) {
   nBaseLevels <- unlist( lapply( mfNoData, FUN=nBaseLevels) )
   stopifnot( inherits(obj, ".cvamEstimateRequest") )
   stopifnot( inherits(model, ".cvamModel") )
   aV <- attr(obj, "vars")
   fV <- aV[ attr(obj, "fixed") ]
   rV <- setdiff(aV, fV)
   form <- as.formula( paste("~", paste(aV, collapse="+") ), 
      env=environment(model) )
   cDT <- as.data.frame( xtabs(form, data=mf0), responseName="freq" )
      cDT$freq <- NULL
   nVarTable <- length(cDT)  # how many variables in cDT
   dimTable <- nBaseLevels[ names(cDT) ]   # dimensions of cDT
   whichVarTable <- match( names(cDT), names(mf0) ) # positions in mf0
   fixedVarTable <- aV %in% fV
   formulaStr <- attr(obj, "formulaStr")
   cDT <- structure( cDT,
      formula = obj,
      formulaStr = formulaStr,
      nVarTable = length(cDT),
      dimTable = dimTable,
      whichVarTable = whichVarTable,
      fixed = structure( fixedVarTable, names=names(cDT) ),
      estimateType = 1L, # 1 = probs
      varConditionedOn = fV,
      varShown = rV,
      varMarginalizedOver = setdiff( attr(model, "vars"), aV),
      class = c( ".cvamEstimate", class(cDT) ) )
   cDT$prob <- 0
   cDT
}

.formatEstimate <- function( obj, confidence, probRound, probDigits,
      fromSeries = FALSE, fromMCMC = FALSE, nDraws = 0L, meanSeries=TRUE ) {
   stopifnot( inherits(obj, ".cvamEstimate") )
   formulaStr <- attr(obj, "formulaStr")
   vco <- attr(obj, "varConditionedOn" )
   vs <- attr(obj, "varShown" )
   SE <- obj$SE
   prob.lower <- obj$prob.lower
   prob.upper <- obj$prob.upper
   obj <- obj[ c(vco, vs, "prob") ]
   if( ( ! is.null(SE) ) && ( ! fromSeries )  ){
      eta <- etaSE <- etaLow <-  etaHigh <- halfWidth <- 
        prob.lower <- prob.upper <- rep( numeric(), length(SE) )
      p <- obj$prob
      w <- ( p > 0 ) & ( p < 1 )
      eta[w] <- log( p[w] / ( 1 - p[w] ) )
      etaSE[w] <- SE[w] / ( p[w] * ( 1 - p[w] ) )
      halfWidth[w] <- etaSE[w] * qnorm( 1 - (1 - confidence)/2 )
      etaLow[w] <- eta[w] - halfWidth[w] 
      etaHigh[w]<- eta[w] + halfWidth[w]
      prob.lower[w] <- 1 / ( 1 + exp( - etaLow[w] ) )
      prob.upper[w] <- 1 / ( 1 + exp( - etaHigh[w] ) )
      w <- ( SE == 0 )
      SE[w] <- prob.lower[w] <- prob.upper[w] <- NA
      obj$SE <- SE
      obj$prob.lower <- prob.lower
      obj$prob.upper <- prob.upper
   }
   obj$SE <- SE
   obj$prob.lower <- prob.lower
   obj$prob.upper <- prob.upper
   if( probRound ) {
      obj$prob <- round(obj$prob, digits=probDigits)
      if( !is.null(obj$SE) ) 
         obj$SE <- round(obj$SE, digits=probDigits)
      if( !is.null(obj$prob.lower) ) 
         obj$prob.lower <- round(obj$prob.lower, digits=probDigits)
      if( !is.null(obj$prob.upper) ) 
         obj$prob.upper <- round(obj$prob.upper, digits=probDigits)
   }
   attr(obj, "nVarTable") <- attr(obj, "dimTable") <- 
      attr(obj, "whichVarTable") <- attr(obj, "fixed") <- 
      attr(obj, "estimateType") <- attr(obj, "varConditionedOn") <- 
      attr(obj, "varShown") <- attr(obj, "varMarginalizedOver") <- NULL
   class(obj) <- c("cvamEstimate", "data.frame")
   attr(obj, "formulaStr") <- formulaStr
   attr(obj, "fromSeries") <- fromSeries
   header <- NULL 
   if( fromSeries ) {
      header <- gettextf(
         "Direct estimates and SEs from %i samples", nDraws )
   } else {
      if( fromMCMC ) {
         if( is.null(SE) ) {
            header <- if( meanSeries ) gettextf(
            "Estimates from %i samples", nDraws ) else
            gettext( "Parameters based on final simulated draw " )
         } else {
            header <- gettextf(
            "Estimates and SEs from %i samples, linearized",
            nDraws )
         }
      } else {
         header <- if( is.null(SE) ) gettext( "Estimates from EM" ) else
            gettext( "Estimates and SE's from EM, linearized" )
      }
   }
   attr(obj, "header") <- header
   obj
}

.cvamEstimateList <- function( obj ) {
   if( inherits(obj, "cvamEstimate") ) return(obj)
   stopifnot( is.list(obj) )
   stopifnot( inherits(obj[[1L]], "cvamEstimate") )
   attr( obj, "header" ) <- attr( obj[[1L]], "header" )
   class( obj ) <- c("cvamEstimateList", class(obj) )
   obj
}

print.cvamEstimateList <- function( x, showHeader=TRUE, ... ){
   stopifnot( inherits(x, "cvamEstimateList") )
   if( showHeader ) cat( attr(x, "header" ), sep="\n" )
   for(i in seq_along(x) ) print( x[[i]], showHeader=FALSE )
   invisible(x)
}

.cvamEstimateSummaries <- function( est ){
   # est is a list of cvamEstimate objects
   ncells <- lapply( est, attr, "dimTable" )
   ncells <- unlist( lapply( ncells, prod ) )
   nEstimates <- length(est)
   estimateType <- unlist( lapply( est, attr, "estimateType") )
   nVarEstimate <- unlist( lapply( est, attr, "nVarTable" ) )
   nVarEstimateTot <- sum( nVarEstimate )
   dimEstimate <- unlist( lapply( est, attr, "dimTable" ) )
   whichVarEstimate <- unlist( lapply( est, attr, "whichVarTable") )
   fixedVarEstimate <- unlist( lapply( est, attr, "fixed") ) # logical
   fixedVarEstimate <- as.integer(fixedVarEstimate) + 1L
   nPackedEstimates <- sum(ncells)
   finPackedEstimates <- cumsum(ncells)
   stPackedEstimates <- c(0L, finPackedEstimates[-length(ncells)] ) + 1L
   estimateInfo <- if( nEstimates > 0L )
      cbind(
         estimateType = estimateType,              # col 1
         nVarEstimate = nVarEstimate,              # col 2
         stPackedEstimates = stPackedEstimates,    # col 3
         finPackedEstimates = finPackedEstimates   # col 4
      ) else matrix( integer(), 0L, 4L )
   storage.mode( estimateInfo ) <- "integer"
   estimateVarInfo <- if( nEstimates > 0L )
      cbind(
         dimEstimate = dimEstimate,             # col 1
         whichVarEstimate = whichVarEstimate,   # col 2
         fixedVarEstimate = fixedVarEstimate    # col 3
      ) else matrix( integer(), 0L, 3L )
   storage.mode( estimateVarInfo ) <- "integer"
   list(
      nEstimates = as.integer(nEstimates),
      nVarEstimateTot = as.integer(nVarEstimateTot),
      nPackedEstimates = as.integer(nPackedEstimates),
      estimateInfo = estimateInfo,
      estimateVarInfo = estimateVarInfo )
}

cvamEstimate <- function( estimate, obj, meanSeries = TRUE,
   confidence = obj$control$confidence,
   probRound = obj$control$probRound,
   probDigits = obj$control$probDigits, ...) {
   if( is.list(estimate) ){
      ios <- unlist( lapply( estimate, .isOneSided ) )
      if( any( !ios ) ) stop( gettext(
         "Some elements of 'estimate' are not one-sided formulas"),
         domain = NA )
   } else {
      if( ! .isOneSided(estimate) ) stop( gettext(
         "'estimate' is not a one-sided formula"), domain = NA )
      estimate <- list(estimate)
   }
   stopifnot( inherits(obj, "cvam") )
   meanSeries <- as.logical(meanSeries)[1L]
   #--------------------------------------------
   # process estimate formulas wrt model and data
   est <- lapply( estimate, .cvamEstimateRequest, obj$model )
   mfNoData <- obj$mfSeen[integer(),,drop=FALSE]
   mfNoData$freq <- NULL
   mf0 <- data.frame( lapply( mfNoData, FUN=dropCoarseLevels ) )
   est <- lapply( est, .cvamEstimate, obj$model, mfNoData, mf0 )
   estSumm <- .cvamEstimateSummaries(est)
   #--------------------------------------------
   nPackedSEs <- if( attr(obj$model, "saturated") ) 0L else
      estSumm$nPackedEstimates
   dimVecEst <- c(
      nEstimates = estSumm$nEstimates,                       # 1
      nVarEstimateTot = estSumm$nVarEstimateTot,             # 2
      nPackedEstimates = estSumm$nPackedEstimates,           # 3
      nPackedSEs = nPackedSEs )                              # 4
   storage.mode( dimVecEst ) <- "integer"
   #--------------------------------------------
   # create a matrix for holding message codes
   msg.len.max <- 40L
   msg.codes <- matrix( 0L, msg.len.max, 17L )
   #--------------------------------------------
   # create dummy input data and frequencies with zero rows
   inputData <- obj$inputData[integer(),,drop=FALSE]
   inputDataFreqInt <- obj$inputDataFreqInt[integer()]
   dimVec <- obj$dimVec
   dimVec[1L] <- 0L
   #--------------------------------------------
   # prepare phat and vhatBeta, and skipSEsInt
   if( obj$method == "EM" ) {
      prob <- obj$prob
      vhatBeta <- obj$vhatBeta
      skipSEsInt <- 0L
   } else {
      if( meanSeries ) {
         prob <- obj$probMean
         vhatBeta <- obj$betaCovMat
         skipSEsInt <- 0L
      } else {
         prob <- obj$prob
         vhatBeta <- obj$betaCovMat   # will be ignored
         skipSEsInt <- 1L
      }
   }
   if( any( ! is.finite(vhatBeta) ) ) vhatBeta[] <- 0
   storage.mode(vhatBeta) <- "double"
   tmp <- .Fortran( "cvam_estimate_em",
      # inputs
      modelTypeInt = obj$modelTypeInt,         
      methodInt = obj$methodInt,
      dimVec = dimVec,
      inputData = inputData,
      inputDataFreqInt = inputDataFreqInt,
      nLevelsMatrix = obj$nLevelsMatrix,
      packedMap = obj$packedMap,
      modelMatrix = obj$modelMatrix,
      offset = obj$offset,
      strZeroInt = obj$strZeroInt,
      prob = prob,
      beta = obj$beta,
      vhatBeta = vhatBeta,
      dimVecEst = dimVecEst,
      estimateInfo = estSumm$estimateInfo,
      estimateVarInfo = estSumm$estimateVarInfo,
      skipSEsInt = skipSEsInt,
      # outputs
      packedEstimates = numeric(estSumm$nPackedEstimates),
      packedSEs = numeric(nPackedSEs),
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
   }
   else{
      msg <- paste0( msg.lines, collapse="\n" )
   }
   msg <- paste( msg, "\n", sep="")
   if( msg!= "OK\n" ) cat( paste("Note: ", msg, sep="") )
   #--------------------------------------------
   if( tmp$status != 0 ) stop( gettext( 
      "Procedure aborted" ), domain = NA )
   #--------------------------------------------
   for( i in seq_along(est) ) {
      st <- estSumm$estimateInfo[i,3L]
      fin <- estSumm$estimateInfo[i,4L]
      est[[i]]$prob <- tmp$packedEstimates[st:fin]
      if( ( length(tmp$packedSEs) > 0L ) & ( skipSEsInt == 0L ) )
         est[[i]]$SE <- tmp$packedSEs[st:fin]
   }
   nDraws <- if( obj$method == "MCMC" | obj$method == "approxBayes") 
      obj$nSampleActual else 0L
   est <- lapply( est, .formatEstimate, confidence,
      probRound, probDigits, fromSeries = FALSE, 
      fromMCMC = ( obj$method == "MCMC" | obj$method == "approxBayes" ),
      nDraws = nDraws, meanSeries = meanSeries )
   est <- if( length(est) == 1L ) est[[1L]] else est
   est <- .cvamEstimateList(est)
   return(est)
}

print.cvamEstimate <- function(x, showHeader = TRUE, ...) {
   stopifnot( inherits(x, "cvamEstimate") )
   if( showHeader ) cat( attr(x, "header"), sep="\n" )
   cat( attr(x, "formulaStr"), sep="\n" )
   print.data.frame(x, ...)
}

anova.cvam <- function( object, ..., method=c("lrt", "logP", "AIC", "BIC"),
   pval=FALSE, pvalDigits=4L, showRank=NULL ) {
   method <- match.arg(method)
   dotargs <- list(...)
   named <- if (is.null(names(dotargs))) 
        rep_len(FALSE, length(dotargs))
   else (names(dotargs) != "")
   if (any(named)) warning(
      "the following arguments to 'anova.glm' are invalid and dropped: ", 
      paste(deparse(dotargs[named]), collapse = ", ") )
   dotargs <- dotargs[!named]
   modList <- c( list(object), dotargs )
   if( length(modList) < 2L ) stop( gettext(
      'Need at least two objects of class "cvam" to compare'), domain = NA ) 
   is.cvam <- vapply(modList, function(x) inherits(x, "cvam"), NA)
   if( any( !is.cvam ) ) stop( gettext(
      'Some supplied objects are not of class "cvam"'), domain = NA ) 
   summList <- lapply( modList, summary.cvam )
   is.EM <- unlist( lapply( summList, `[[`, "method" ) ) == "EM"
   anyLatent <- any( unlist( lapply( summList, `[[`, "anyLatent" ) ) )
   if( any( !is.cvam ) ) warning( gettext(
      'Some supplied objects do not have method "EM"'), domain = NA ) 
   nCells <- unlist( lapply( summList, `[[`, "nCells" ) )
   nCellsNonLatent <- unlist( lapply( summList, `[[`, "nCellsNonLatent" ) )
   if( any( nCellsNonLatent != nCellsNonLatent[1L] ) ) warning( gettext(
      "Models do not have the same number of cells"), domain = NA )
   nZero <- unlist( lapply( summList, `[[`, "nZero" ) )
   if( any( nZero != nZero[1L] ) ) warning( gettext(
      "Models do not have the same number of structural zeros"), domain = NA )
   nTotal <- unlist( lapply( summList, `[[`, "totalSampSize" ) )
   if( any( nTotal != nTotal[1L] ) ) warning( gettext(
      "Fitted models are based on different sample sizes"), domain = NA )
   priorSampSize <- unlist( lapply( summList, `[[`, "priorSampSize" ) )
   if( any( priorSampSize != priorSampSize[1L] ) ) warning( gettext(
      "Fitted models have different prior sample sizes"), domain = NA )
   mList <- lapply( summList, `[[`, "model" )
   formulaStr <- unlist( lapply( mList, attr, "formulaStr" ) )
   formulaStr <- paste("Model ", format(1:length(summList)), ": ",
      formulaStr, sep="")
   formulaStr <- paste( formulaStr, collapse="\n" )
   resid.df <- unlist( lapply( summList, `[[`, "degreesOfFreedom" ) )
   if( method %in% c("lrt", "logP") ) {
      meas <- if( method == "lrt" ) 
         unlist( lapply( summList, `[[`, "loglik" ) ) else
         unlist( lapply( summList, `[[`, "logP" ) )
      meas <- -2*meas
      result <- data.frame( resid.df, meas )
      names(result)[2L] <- if( method == "lrt" ) "-2*loglik" else "-2*logP" 
      result$df <- c( NA, resid.df[-length(resid.df)]) - resid.df
      result$change <- c( NA, meas[-length(meas)] ) - meas
      pvalDigits <- as.integer(pvalDigits)[1L]
      if( pval ) result$pval <- 
         round( 1 - pchisq( result$change, result$df ), pvalDigits )
   } else {
      meas <- -2 * unlist( lapply( summList, `[[`, "loglik" ) )
      result <- data.frame( resid.df, meas )
      names(result)[2] <- "-2*loglik" 
      nParams <- unlist( lapply( summList, `[[`, "nParams" ) )
      IC <- if( method == "AIC" ) meas + 2*nParams else
         meas + log(nTotal) * nParams
      if( method == "AIC" ) meas <- result$AIC <- IC else
         meas <- result$BIC <- IC
   }
   showRank <- if( is.null(showRank) ) method %in% c("AIC", "BIC") else
      as.logical(showRank)[1L] 
   if( showRank ) result$rank <- rank(meas)
   structure( result,
      heading = formulaStr,
      class = c("anova", "data.frame") )
}

cvamPredict <- function( form, obj, data, freq, meanSeries = TRUE, sep="." ){
   if( !.isOneSided(form) ) stop( gettext(
      "'form' is not a one-sided formula"), domain = NA)
   stopifnot( inherits(obj, "cvam") )
   sep <- as.character(sep)[1L]
   #----------------------------------
   pred <- .cvamPredictRequest(form, obj$model )
   #----------------------------------
   # get the prediction data frame and frequencies
   mc <- match.call( expand.dots=FALSE )
   mc[[1]] <- quote(stats::model.frame)
   m <- match( c("form", "data", "freq"), names(mc), nomatch=0L )
   mc <- mc[ c(1L,m) ]
   names(mc)[2] <- "formula"
   mc[[2]] <- obj$model
   mc$na.action <- as.name("na.pass")
   mc$drop.unused.levels <- FALSE
   predF <- eval( mc, parent.frame() )
   rN <- rownames(predF)
   if( is.null( predF$`(freq)` ) ) {
      freq <- rep(1L, NROW(predF))
      freqSupplied <- FALSE
   } else {
      freq <- predF$`(freq)`
      predF$`(freq)` <- NULL 
      freqSupplied <- TRUE
   }
   freqInt <- as.integer(freq)
   if( any( freq != freqInt ) ) warning( gettext(
      "Some frequencies changed when integerized" ), domain=NA )
   badV <- names(predF)[ ! unlist( lapply( predF, is.factor ) ) ]
   if( length(badV) > 0L ) stop( gettextf(
      "Variable '%s' is not a factor", badV[1L] ), domain=NA )
   if( any(is.na(freq)) ) stop( gettext(
      "Missing values in 'freq' not allowed"), domain = NA )
   predF <- data.frame( lapply( predF, FUN=coarsened, warnIfCoarsened=FALSE) )
#   predF <- data.frame( lapply( predF, FUN=sticky::sticky ) )
   #----------------------------------
   # check for consistency with obj
   predFNoData <- predF[integer(),,drop=FALSE]
   predF0 <- data.frame( lapply( predFNoData, FUN=dropCoarseLevels ) )
   mfNoData <- obj$mfSeen[integer(),,drop=FALSE]
   mfNoData$freq <- NULL
   mf0 <- data.frame( lapply( mfNoData, FUN=dropCoarseLevels ) )
   if( !all.equal( names(predF0), names(mf0) ) ) stop( gettext(
      "Wrong variable names in prediction frame"), domain = NA)
   for(i in seq_along(predF0) ) {
      if( !identical( attributes(predF0[[i]]), attributes(mf0[[i]]) ) )
         stop( gettextf(
         "Attributes of variable '%s' in prediction frame", names(predF)[i] ),
         gettextf("do not match those of variable '%s' in 'obj'",
            names(predF)[i] ), domain = NA)
   }
   #--------------------------------------------
   # check for missing or coarsened values in fixed variables 
   if( any( attr(obj$model, "fixed") ) ) {
      mfFixed <- predF[ attr(obj$model, "fixed") ]
      mfFixed <- lapply( mfFixed, dropCoarseLevels )
      mfFixed <- lapply( mfFixed, is.na )
      mfFixed <- unlist( lapply( mfFixed, any ) )
      if( any(mfFixed) ) { 
         fV <- attr(obj$model, "vars")[ attr(obj$model, "fixed") ]
         fVbad <- fV[mfFixed]
         stop( gettextf( 
            "Non-modeled variable '%s' contains missing or coarsened values",
            fVbad[1L] ), domain = NA ) 
      }
   }
   #--------------------------------------------
   # process predict formula wrt model and data
   pred <- .cvamPredict( pred, obj$model, predF, predFNoData, predF0, sep )
   predSumm <- .cvamPredictSummaries(pred)
   #--------------------------------------------
   # prepare for Fortran
   predData <- data.matrix(predF)
   storage.mode(predData) <- "integer"
   dimVec <- obj$dimVec
   dimVec[1L] <- NROW(predData)
   predMat <- as.matrix(pred)
   storage.mode(predMat) <- "double"
   pI <- predSumm$predictInfo
   dimVecPred <- c(
      nRowPredData = NROW(predData),       # 1
      nCells = pI[["nCells"]],             # 2
      nVarPredict = pI[["nVarPredict"]] )    # 3
   vhatBeta <- obj$vhatBeta
   if( any( ! is.finite(vhatBeta) ) ) vhatBeta[] <- 0
   storage.mode(vhatBeta) <- "double"
   #--------------------------------------------
   # which probabilities to use
   if( obj$method == "EM" ) {
      prob <- obj$prob
   } else {
      prob <- if( meanSeries ) obj$probMean else obj$prob
   }
   #--------------------------------------------
   # create a matrix for holding message codes
   msg.len.max <- 40L
   msg.codes <- matrix( 0L, msg.len.max, 17L )
   #--------------------------------------------
   tmp <- .Fortran( "cvam_predict_em",
      # inputs
      modelTypeInt = obj$modelTypeInt,
      methodInt = obj$methodInt,
      dimVec = dimVec,
      predData = predData,
      predDataFreqInt = freqInt,
      nLevelsMatrix = obj$nLevelsMatrix,
      packedMap = obj$packedMap,
      modelMatrix = obj$modelMatrix,
      offset = obj$offset,
      strZeroInt = obj$strZeroInt,
      prob = prob,
      beta = obj$beta,
      vhatBeta = vhatBeta,
      dimVecPred = dimVecPred,
      predictVarInfo = predSumm$predictVarInfo,
      # outputs
      predMat = predMat,
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
   }
   else{
      msg <- paste0( msg.lines, collapse="\n" )
   }
   msg <- paste( msg, "\n", sep="")
   if( msg!= "OK\n" ) cat( paste("Note: ", msg, sep="") )
   #--------------------------------------------
   if( tmp$status != 0 ) stop( gettext( 
      "Procedure aborted" ), domain = NA )
   #--------------------------------------------
   # process predicts         
   result <- structure( as.data.frame( tmp$predMat ),
      colFrame=attr(pred,"colFrame"),
      class = c("cvamPredict", "data.frame") )
   rownames(result) <- rN
   return(result)
}


.cvamPredictRequest <- function( obj, model ){
   aN <- setdiff( all.names(obj), c("~", "(", "+") )
   allVars <- aV <- all.vars(obj)
   badN <- setdiff( aN, aV )
   if( length(badN) > 0L ) stop( gettextf(
      "In formula '%s': symbol '%s' not allowed", 
         deparse(obj), badN[1L] ), domain = NA )
   notInModel <- setdiff( allVars, attr(model, "vars") )
   if( length(notInModel) > 0L ) stop( gettextf(
      "In formula '%s': variable '%s' does not appear in the model",
      deparse(obj),  notInModel[1L] ), domain = NA )
   formulaStr <- deparse(obj)
   formulaStr <- paste("~ ", substring(formulaStr, 2L), sep="" )
   structure( obj,
      vars = allVars,
      formulaStr = formulaStr,
      .Environment = attr(model, ".Environment"),
      class = c( ".cvamPredictRequest", class(obj) ) )
}

.cvamPredict <- function( obj, model, predF, predFNoData, predF0, sep ) {
   nBaseLevels <- unlist( lapply( predFNoData, FUN=nBaseLevels) )
   stopifnot( inherits(obj, ".cvamPredictRequest") )
   stopifnot( inherits(model, ".cvamModel") )
   aV <- attr(obj, "vars")
   levs <- if( length(aV) == 1L ) levels( predF0[[aV]] ) else
      levels( interaction( predF0[aV], sep=sep ) )
   pred <- matrix( numeric(1L), NROW(predF), length(levs) )
   rownames(pred) <- rownames(predF)
   colnames(pred) <- levs
   colFrame <- as.data.frame( table( predF0[aV] ), responseName = "freq" )
   colFrame$freq <- NULL
   structure( as.data.frame(pred),
      formula = obj,
      formulaStr = attr(obj, "formulaStr"),
      nVarTable = length(aV),
      dimTable = nBaseLevels[aV],
      whichVarTable = match(aV, names(predF0) ),
      colFrame = colFrame,
      class = c( ".cvamPredict", "data.frame" ) )
}

.cvamPredictSummaries <- function( pred ) {
   nCells <- prod( attr(pred, "dimTable") )
   nVarPredict <- attr(pred, "nVarTable")
   dimPredict <- attr(pred, "dimTable")
   whichVarPredict <- attr(pred, "whichVarTable")
   predictInfo <- c( nCells=nCells, nVarPredict=nVarPredict )
   storage.mode(predictInfo) <- "integer"
   predictVarInfo <- cbind( dimPredict = dimPredict,        # col 1
                       whichVarPredict = whichVarPredict )  # col 2
   storage.mode(predictVarInfo) <- "integer"
   list( predictInfo = predictInfo,
      predictVarInfo = predictVarInfo )
}

cvamImpute <- function( obj, data, freq, meanSeries = FALSE, synthetic=FALSE ){
   stopifnot( inherits(obj, "cvam") )
   #----------------------------------
   # get the imputation data frame and frequencies
   mc <- match.call( expand.dots=FALSE )
   mc[[1]] <- quote(stats::model.frame)
   m <- match( c("obj", "data", "freq"), names(mc), nomatch=0L )
   mc <- mc[ c(1L,m) ]
   names(mc)[2] <- "formula"
   mc[[2]] <- obj$model
   mc$na.action <- as.name("na.pass")
   mc$drop.unused.levels <- FALSE
   impF <- eval( mc, parent.frame() )
   rN <- rownames(impF)
   if( is.null( impF$`(freq)` ) ) {
      freq <- rep(1L, NROW(impF))
      freqSupplied <- FALSE
   } else {
      freq <- impF$`(freq)`
      impF$`(freq)` <- NULL 
      freqSupplied <- TRUE
   }
   freqInt <- as.integer(freq)
   if( any( freq != freqInt ) ) warning( gettext(
      "Some frequencies changed when integerized" ), domain=NA )
   badV <- names(impF)[ ! unlist( lapply( impF, is.factor ) ) ]
   if( length(badV) > 0L ) stop( gettextf(
      "Variable '%s' is not a factor", badV[1L] ), domain=NA )
   if( any(is.na(freq)) ) stop( gettext(
      "Missing values in 'freq' not allowed"), domain = NA )
   impF <- data.frame( lapply( impF, FUN=coarsened, warnIfCoarsened=FALSE) )
#   impF <- data.frame( lapply( impF, FUN=sticky::sticky ) )
   #----------------------------------
   # check for consistency with obj
   impFNoData <- impF[integer(),,drop=FALSE]
   impF0 <- data.frame( lapply( impFNoData, FUN=dropCoarseLevels ) )
   mfNoData <- obj$mfSeen[integer(),,drop=FALSE]
   mfNoData$freq <- NULL
   mf0 <- data.frame( lapply( mfNoData, FUN=dropCoarseLevels ) )
   if( !all.equal( names(impF0), names(mf0) ) ) stop( gettext(
      "Wrong variable names in imputation frame"), domain = NA)
   for(i in seq_along(impF0) ) {
      if( !identical( attributes(impF0[[i]]), attributes(mf0[[i]]) ) )
         stop( gettextf(
         "Attributes of variable '%s' in imputation frame", names(impF)[i] ),
         gettextf("do not match those of variable '%s' in 'obj'",
            names(impF)[i] ), domain = NA)
   }
   #--------------------------------------------
   # check for missing or coarsened values in fixed variables 
   if( any( attr(obj$model, "fixed") ) ) {
      mfFixed <- impF[ attr(obj$model, "fixed") ]
      mfFixed <- lapply( mfFixed, dropCoarseLevels )
      mfFixed <- lapply( mfFixed, is.na )
      mfFixed <- unlist( lapply( mfFixed, any ) )
      if( any(mfFixed) ) { 
         fV <- attr(obj$model, "vars")[ attr(obj$model, "fixed") ]
         fVbad <- fV[mfFixed]
         stop( gettextf( 
            "Non-modeled variable '%s' contains missing or coarsened values",
            fVbad[1L] ), domain = NA ) 
      }
   }
   if( freqSupplied ) {
      #--------------------------------------------
      # aggregate the input data
      form <-  as.formula( paste( "freq ~", 
         paste( attr(obj$model, "vars"), collapse="+" ) ), 
         env = environment(obj$model) )
      impF$freq <- freqInt
      impF <- aggregate(form, FUN=sum, data=impF)
      freqInt <- impF$freq
      impF$freq <- NULL
      #--------------------------------------------
      # create result with correct size, shape and attributes
      nCells <- prod( unlist( lapply( impF0, nlevels ) ) )
      result <- impF0[ 1:nCells,,drop=FALSE]
      rownames(result) <- NULL
      resultMat <- data.matrix(result)
      resultMat[] <- 0L
      storage.mode(resultMat) <- "integer"
      resultFreqInt <- integer( NROW(resultMat) )
      #--------------------------------------------
      # prepare for Fortran call
      inputData <- data.matrix(impF)
      storage.mode(inputData) <- "integer"
      inputDataFreqInt <- freqInt
      storage.mode(inputDataFreqInt) <- "integer"
      dimVec <- obj$dimVec
      dimVec[1L] <- NROW(inputData)
      vhatBeta <- obj$vhatBeta
      if( any( ! is.finite(vhatBeta) ) ) vhatBeta[] <- 0
      storage.mode(vhatBeta) <- "double"
      #--------------------------------------------
      # which probabilities to use
      if( obj$method == "EM" ) {
         prob <- obj$prob
      } else {
         prob <- if( meanSeries ) obj$probMean else obj$prob
      }
      #--------------------------------------------
      # create a matrix for holding message codes
      msg.len.max <- 40L
      msg.codes <- matrix( 0L, msg.len.max, 17L )
      #--------------------------------------------
      tmp <- .Fortran( "cvam_impute_freq",
         # inputs
         modelTypeInt = obj$modelTypeInt,
         methodInt = obj$methodInt,
         dimVec = dimVec,
         inputData = inputData,
         inputDataFreqInt = inputDataFreqInt,
         nLevelsMatrix = obj$nLevelsMatrix,
         packedMap = obj$packedMap,
         modelMatrix = obj$modelMatrix,
         offset = obj$offset,
         strZeroInt = obj$strZeroInt,
         prob = prob,
         beta = obj$beta,
         vhatBeta = vhatBeta,
         syntheticInt = as.integer( as.logical(synthetic)[1L] ),
         # outputs
         resultMat = resultMat,
         resultFreqInt = resultFreqInt,
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
      }
      else{
         msg <- paste0( msg.lines, collapse="\n" )
      }
      msg <- paste( msg, "\n", sep="")
      if( msg!= "OK\n" ) cat( paste("Note: ", msg, sep="") )
      #--------------------------------------------
      if( tmp$status != 0 ) stop( gettext( 
         "Procedure aborted" ), domain = NA )
      #--------------------------------------------
      for( i in 1:length(result) ) {
         levs <- levels( result[[i]] )
         result[[i]] <- factor( tmp$resultMat[,i],
            ordered = is.ordered( result[[i]] )  )
         levels( result[[i]] ) <- levs
      }
      result$freq <- tmp$resultFreqInt
      return(result)
   } else {
      #--------------------------------------------
      # keep impF as microdata
      # create result with correct size, shape and attributes
      result <- data.frame( lapply(impF, FUN=dropCoarseLevels) )
      rownames(result) <- rN
      #--------------------------------------------
      # prepare for Fortran call
      inputData <- data.matrix(impF)
      storage.mode(inputData) <- "integer"
      inputDataFreqInt <- freqInt
      storage.mode(inputDataFreqInt) <- "integer"
      dimVec <- obj$dimVec
      dimVec[1L] <- NROW(inputData)
      vhatBeta <- obj$vhatBeta
      if( any( ! is.finite(vhatBeta) ) ) vhatBeta[] <- 0
      storage.mode(vhatBeta) <- "double"
      resultMat <- inputData
      #--------------------------------------------
      # which probabilities to use
      if( obj$method == "EM" ) {
         prob <- obj$prob
      } else {
         prob <- if( meanSeries ) obj$probMean else obj$prob
      }
      #--------------------------------------------
      # create a matrix for holding message codes
      msg.len.max <- 40L
      msg.codes <- matrix( 0L, msg.len.max, 17L )
      #--------------------------------------------
      tmp <- .Fortran( "cvam_impute_microdata",
         # inputs
         modelTypeInt = obj$modelTypeInt,
         methodInt = obj$methodInt,
         dimVec = dimVec,
         inputData = inputData,
         inputDataFreqInt = inputDataFreqInt,
         nLevelsMatrix = obj$nLevelsMatrix,
         packedMap = obj$packedMap,
         modelMatrix = obj$modelMatrix,
         offset = obj$offset,
         strZeroInt = obj$strZeroInt,
         prob = prob,
         beta = obj$beta,
         vhatBeta = vhatBeta,
         syntheticInt = as.integer( as.logical(synthetic)[1L] ),
         # outputs
         resultMat = resultMat,
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
      }
      else{
         msg <- paste0( msg.lines, collapse="\n" )
      }
      msg <- paste( msg, "\n", sep="")
      if( msg!= "OK\n" ) cat( paste("Note: ", msg, sep="") )
      #--------------------------------------------
      if( tmp$status != 0 ) stop( gettext( 
         "Procedure aborted" ), domain = NA )
      #--------------------------------------------
      for( i in 1:length(result) ) {
         levs <- levels( result[[i]] )
         result[[i]] <- factor( tmp$resultMat[,i],
            ordered = is.ordered( result[[i]] )  )
         levels( result[[i]] ) <- levs
      }
      return(result)
   }
}

cvamLik <- function( form, obj, data, meanSeries = TRUE ) {
   if( !.isOneSided(form) ) stop( gettext(
      "'form' is not a one-sided formula"), domain = NA)
   stopifnot( inherits(obj, "cvam") )
   #--------------------------------------------
   # process formula wrt model and data
   lik0 <- .cvamLikRequest(form, obj$model)
   mfNoData <- obj$mfSeen[integer(),,drop=FALSE]
   mfNoData$freq <- NULL
   mf0 <- data.frame( lapply( mfNoData, FUN=dropCoarseLevels ) )
   lik <- .cvamLik( lik0, obj$model, mfNoData, mf0 )
   #--------------------------------------------
   # create data frame
   mc <- match.call( expand.dots=FALSE )
   mc[[1]] <- quote(stats::model.frame)
   m <- match( c("form", "data"), names(mc), nomatch=0L )
   mc <- mc[ c(1L,m) ]
   names(mc)[2] <- "formula"
   mc[[2]] <- obj$model
   mc$na.action <- as.name("na.pass")
   mc$drop.unused.levels <- FALSE
   likF <- eval( mc, parent.frame() )
   rN <- rownames(likF)
   badV <- names(likF)[ ! unlist( lapply( likF, is.factor ) ) ]
   if( length(badV) > 0L ) stop( gettextf(
      "Variable '%s' is not a factor", badV[1L] ), domain=NA )
   likF <- data.frame( lapply( likF, FUN=coarsened, warnIfCoarsened=FALSE) )
#   likF <- data.frame( lapply( likF, FUN=sticky::sticky ) )
   #--------------------------------------------
   # check for consistency with obj
   likFNoData <- likF[integer(),,drop=FALSE]
   likF0 <- data.frame( lapply( likFNoData, FUN=dropCoarseLevels ) )
   mfNoData <- obj$mfSeen[integer(),,drop=FALSE]
   mfNoData$freq <- NULL
   mf0 <- data.frame( lapply( mfNoData, FUN=dropCoarseLevels ) )
   if( !all.equal( names(likF0), names(mf0) ) ) stop( gettext(
      "Wrong variable names in prediction frame"), domain = NA)
   for(i in seq_along(likF0) ) {
      if( !identical( attributes(likF0[[i]]), attributes(mf0[[i]]) ) )
         stop( gettextf(
         "Attributes of variable '%s' in data frame", names(likF)[i] ),
         gettextf("do not match those of variable '%s' in 'obj'",
            names(likF)[i] ), domain = NA)
   }
   #--------------------------------------------
   # check for missing or coarsened values in variables fixed by model 
   if( any( attr(obj$model, "fixed") ) ) {
      mfFixed <- likF[ attr(obj$model, "fixed") ]
      mfFixed <- lapply( mfFixed, dropCoarseLevels )
      mfFixed <- lapply( mfFixed, is.na )
      mfFixed <- unlist( lapply( mfFixed, any ) )
      if( any(mfFixed) ) { 
         fV <- attr(obj$model, "vars")[ attr(obj$model, "fixed") ]
         fVbad <- fV[mfFixed]
         stop( gettextf( 
            "Non-modeled variable '%s' contains missing or coarsened values",
            fVbad[1L] ), domain = NA ) 
      }
   }
   #--------------------------------------------
   # check for coarsened values in variables fixed in lik formula 
   if( any( attr(lik0, "fixed") ) ) {
      fV <- attr(lik0, "vars")[ attr(lik0, "fixed") ]
      mfFixed <- likF[ fV ]
      mfFixed <- lapply( mfFixed, dropCoarseLevels )
      mfFixed <- lapply( mfFixed, is.na )
      mfFixed <- unlist( lapply( mfFixed, any ) )
      if( any(mfFixed) ) { 
         fVbad <- fV[mfFixed]
         stop( gettextf( 
            "Variable '%s', which is fixed in 'form', ", fVbad[1L] ),
            gettext( "has missing or coarsened values" ), domain = NA )
      }
   }
   #--------------------------------------------
   # prepare for Fortran
   likSumm <- .cvamEstimateSummaries( list(lik) )
   likVarInfo <- likSumm$estimateVarInf
   likInfo <- c( nCells = prod( likVarInfo[,1] ), # 1
      nVar = NROW(likVarInfo) )                   # 2
   inputData <- data.matrix(likF)
   storage.mode(inputData) <- "integer"
   inputDataFreqInt <- integer( NROW(inputData) ) # set to zero
   dimVec <- obj$dimVec
   dimVec[1L] <- NROW(inputData)
   likValues <- numeric( NROW(inputData) )
   storage.mode(likValues) <- "double"
   dimVecLik <- likInfo
   storage.mode(dimVecLik) <- "integer"
   vhatBeta <- obj$vhatBeta
   if( any( ! is.finite(vhatBeta) ) ) vhatBeta[] <- 0
   storage.mode(vhatBeta) <- "double"
   #--------------------------------------------
   # which probabilities to use
   if( obj$method == "EM" ) {
      prob <- obj$prob
   } else {
      prob <- if( meanSeries ) obj$probMean else obj$prob
   }
   #--------------------------------------------
   # create a matrix for holding message codes
   msg.len.max <- 40L
   msg.codes <- matrix( 0L, msg.len.max, 17L )
   #--------------------------------------------
   tmp <- .Fortran( "cvam_lik",
      # inputs
      modelTypeInt = obj$modelTypeInt,
      methodInt = obj$methodInt,
      dimVec = dimVec,
      inputData = inputData,
      inputDataFreqInt = inputDataFreqInt,
      nLevelsMatrix = obj$nLevelsMatrix,
      packedMap = obj$packedMap,
      modelMatrix = obj$modelMatrix,
      offset = obj$offset,
      strZeroInt = obj$strZeroInt,
      prob = prob,
      beta = obj$beta,
      vhatBeta = vhatBeta,
      dimVecLik = dimVecLik,
      likVarInfo = likVarInfo,
      # outputs
      likValues = likValues,
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
   }
   else{
      msg <- paste0( msg.lines, collapse="\n" )
   }
   msg <- paste( msg, "\n", sep="")
   if( msg!= "OK\n" ) cat( paste("Note: ", msg, sep="") )
   #--------------------------------------------
   if( tmp$status != 0 ) stop( gettext( 
      "Procedure aborted" ), domain = NA )
   #--------------------------------------------
   likF$likVal <- tmp$likValues
   likF
}

.cvamLikRequest <- function( obj, model ) {
   aN <- setdiff( all.names(obj), c("~", "(", "+", "|") )
   aV <- all.vars(obj)
   badN <- setdiff( aN, aV )
   if( length(badN) > 0L ) stop( gettextf(
      "In formula '%s': symbol '%s' not allowed", 
         deparse(obj), badN[1L] ), domain = NA )
   Form <- Formula::Formula(obj)
   nParts <- length( attr(Form, "rhs") )
   if( nParts > 2L ) stop( gettextf( 
     "Multiple '|' found in '%s'", deparse(obj) ), domain = NA )
   randVars <- all.vars( formula(Form, rhs=1L ) )
   fixedVars <- if( nParts == 2L ) 
      all.vars( formula(Form, rhs=2L) )  else character()
   if( length(randVars) == 0L ) {
      if( nParts == 1L ) stop( gettextf(
         "In formula '%s': no variables present",
         deparse(obj) ), domain = NA ) else stop( gettextf(
      "In formula '%s': no variables present before '|'",
      deparse(obj)), domain = NA )
   }
   inBoth <- intersect( randVars, fixedVars )
   if( length(inBoth) > 0L ) stop( gettextf(
      "In formula '%s': variable '%s' appears on both sides of '|'", 
      deparse(obj), inBoth[1L] ), domain = NA )
   allVars <- c( randVars, fixedVars )
   notInModel <- setdiff( allVars, attr(model, "vars") )
   if( length(notInModel) > 0L ) stop( gettextf(
      "In formula '%s': variable '%s' does not appear in the model",
      deparse(obj),  notInModel[1L] ), domain = NA )
   fInModel <- attr(model, "vars")[ attr(model, "fixed") ]
   notModeled <- intersect( randVars, fInModel )
   if( length(notModeled) > 0L ) stop( gettextf(
      "In formula '%s': variable '%s' is fixed in the model",
      deparse(obj),  notModeled[1L] ), domain = NA )
   missingFixedVars <- setdiff( fInModel, fixedVars )
   if( length(missingFixedVars) > 0L ) stop( gettextf(
      "In formula '%s': variable '%s' does not appear after '|'",
      deparse(obj),  missingFixedVars[1L] ), domain = NA )
   formulaStr <- deparse(obj)
   formulaStr <- paste("~ ", substring(formulaStr, 2L), sep="" )
   structure( obj,
      vars = allVars,
      fixed = allVars %in% fixedVars,
      formulaStr = formulaStr,
      .Environment = attr(model, ".Environment"),
      class = c( ".cvamLikRequest", class(obj) ) )
}

.cvamLik <- function( obj, model, mfNoData, mf0 ) {
   nBaseLevels <- unlist( lapply( mfNoData, FUN=nBaseLevels) )
   stopifnot( inherits(obj, ".cvamLikRequest") )
   stopifnot( inherits(model, ".cvamModel") )
   aV <- attr(obj, "vars")
   fV <- aV[ attr(obj, "fixed") ]
   rV <- setdiff(aV, fV)
   form <- as.formula( paste("~", paste(aV, collapse="+") ), 
      env=environment(model) )
   cDT <- as.data.frame( xtabs(form, data=mf0), responseName="freq" )
      cDT$freq <- NULL
   nVarTable <- length(cDT)  # how many variables in cDT
   dimTable <- nBaseLevels[ names(cDT) ]   # dimensions of cDT
   whichVarTable <- match( names(cDT), names(mf0) ) # positions in mf0
   fixedVarTable <- aV %in% fV
   formulaStr <- attr(obj, "formulaStr")
   cDT <- structure( cDT,
      formula = obj,
      formulaStr = formulaStr,
      nVarTable = length(cDT),
      dimTable = dimTable,
      whichVarTable = whichVarTable,
      fixed = structure( fixedVarTable, names=names(cDT) ),
      varConditionedOn = fV,
      varShown = rV,
      varMarginalizedOver = setdiff( attr(model, "vars"), aV),
      class = c( ".cvamLik", class(cDT) ) )
   cDT$prob <- 0
   cDT
}

get.coef <- function( obj, withSE = FALSE, meanSeries = TRUE,
   msgIfNone = TRUE ) {
   stopifnot( inherits(obj, "cvam") )
   if( attr(obj$model, "saturated") ) {
      if( msgIfNone ) message( gettext(
         "Model is saturated, coefficients are not defined" ), domain = NA )
      return( invisible() )
   }
   if( obj$method == "EM" ) {
      coef <- obj$beta
      if( ! withSE ) return(coef)
      SE <- sqrt( diag(obj$vhatBeta) )
      if( all( obj$vhatBeta == 0 ) ) SE[] <- NA
      tstat <- coef / SE
      pval <- 2 * pnorm( - abs(tstat) )
      cN <- names(coef)
      coef <- data.frame( coef, SE, tstat = round(tstat,2),
         pval = round(pval, 4) )
      row.names(coef) <- cN
      attr( coef, "header" ) <- gettext(
         "Estimates from EM, with Hessian-based SEs" ) 
   } else if( obj$method == "MCMC" ) {
      nDraws <- max( obj$nIterActual -  obj$control$burnMCMC, 0L )
      if( nDraws == 0L && meanSeries ) {
         if( msgIfNone ) message( gettext(
            "No MCMC samples available for estimating coefficients" ),
            domain = NA )
         return( invisible() )
      } else {
         coef <- if( meanSeries ) obj$betaMean else obj$beta
         if( ! withSE ) return(coef)
         if( ! meanSeries ) {
            if( msgIfNone ) message( gettext(
               "Should not request SEs when 'meanSeries' is FALSE" ),
               domain = NA )
            return( invisible() )
         }
         SE <- sqrt( diag( obj$betaCovMat ) )
         tstat <- coef / SE
         pval <- 2 * pnorm( - abs(tstat) )
         cN <- names(coef)
         coef <- data.frame( coef, SE, tstat = round(tstat,2),
            pval = round(pval, 4) )
         row.names(coef) <- cN
         attr( coef, "header" ) <- gettextf(
            "Direct estimates and SE's from %i successive MCMC samples",
            nDraws ) 
      }
   } else if( obj$method == "approxBayes" ) {
      nDraws <- obj$nIterActual
      if( nDraws == 0L && meanSeries ) {
         if( msgIfNone ) message( gettext(
            "No samples available for estimating coefficients" ),
            domain = NA )
         return( invisible() )
      } else {
         coef <- if( meanSeries ) obj$betaMean else obj$beta
         if( ! withSE ) return(coef)
         if( ! meanSeries ) {
            if( msgIfNone ) message( gettext(
               "Should not request SEs when 'meanSeries' is FALSE" ),
               domain = NA )
            return( invisible() )
         }
         SE <- sqrt( diag( obj$betaCovMat ) )
         tstat <- coef / SE
         pval <- 2 * pnorm( - abs(tstat) )
         cN <- names(coef)
         coef <- data.frame( coef, SE, tstat = round(tstat,2),
            pval = round(pval, 4) )
         row.names(coef) <- cN
         attr( coef, "header" ) <- gettextf(
            "Direct estimates and SE's from %i successive samples",
            nDraws ) 
      }
   }
   return(coef)
}

get.estimates <- function( obj, msgIfNone = TRUE ) {
   stopifnot( inherits(obj, "cvam") )
   if( is.null(obj$estimates) & msgIfNone ) {
      message( gettext( "No estimates are present" ), domain = NA )
      return( invisible() )
   }
   obj$estimates
}

get.covMat <- function( obj, msgIfNone = TRUE ) {
   stopifnot( inherits(obj, "cvam") )
   if( attr(obj$model, "saturated") ) {
      if( msgIfNone ) message( gettext(
         "Model is saturated, coefficients are not defined" ), domain = NA )
      return( invisible() )
   }
   if( obj$method == "EM" ) {
      result <- obj$vhatBeta
   } else {
      nDraws <- max( obj$nIterActual -  obj$control$burnMCMC, 0L )
      if( nDraws == 0L ) {
         if( msgIfNone ) message( gettext(
            "No MCMC samples available for estimating covMat" ),
            domain = NA )
         return( invisible() )
      }
      result <- obj$betaCovMat
   }
   result
}


get.loglik <- function( obj, msgIfNone = TRUE ) {
   stopifnot( inherits(obj, "cvam") )
   if( obj$iter < 1L ) {
      if( msgIfNone ) message( gettext(
         "No iterations of EM or MCMC were performed" ),
         domain = NA )
      return( invisible() )
   }
   obj$loglik[ obj$iter ]
}

get.logP <- function( obj, msgIfNone = TRUE ) {
   stopifnot( inherits(obj, "cvam") )
   if( obj$iter < 1L ) {
      if( msgIfNone ) message( gettext(
         "No iterations of EM or MCMC were performed" ),
         domain = NA )
      return( invisible() )
   }
   obj$logP[ obj$iter ]
}

get.mfTrue <- function( obj ){
   stopifnot( inherits(obj, "cvam") )
   obj$mfTrue
}

get.modelMatrix <- function( obj, msgIfNone = TRUE ){
   stopifnot( inherits(obj, "cvam") )
   if( attr(obj$model, "saturated") ) {
      if( msgIfNone ) message( gettextf(
         "Model is saturated, there is no modelMatrix" ), domain = NA )
      return( invisible() )
   }
   obj$modelMatrix
}

get.offset <- function( obj, mfTrue = FALSE, msgIfNone = TRUE ){
   stopifnot( inherits(obj, "cvam") )
   if( attr(obj$model, "saturated") ) {
      if( msgIfNone ) message( gettextf(
         "Model is saturated, there is no offset" ), domain = NA )
      return( invisible() )
   }
   if( mfTrue ) {
      result <- obj$mfTrue
      result$offset <- obj$offset
   } else result <- obj$offset
   result
}

get.strZero <- function( obj, mfTrue = FALSE, msgIfNone = TRUE ){
   if( mfTrue ) {
      result <- obj$mfTrue
      result$strZero <- obj$strZero
   } else result <- obj$strZero
   result
}

get.fitted <- function( obj, type=c("prob", "mean", "logMean"),
   mfTrue = TRUE, meanSeries = TRUE, msgIfNone = TRUE ){
   stopifnot( inherits(obj, "cvam") )
   type <- match.arg(type)
   if( obj$method == "MCMC" && meanSeries ) {
      nDraws <- max( obj$nIterActual -  obj$control$burnMCMC, 0L )
      if( nDraws == 0L && meanSeries ) {
         if( msgIfNone ) message( gettext(
            "No MCMC samples available after burn-in for obtaining fit" ),
            domain = NA )
         return( invisible() )
      }
   }
   if( obj$method == "EM" ) {
      prob <- obj$prob
   } else {
      prob <- if( meanSeries ) obj$probMean else obj$prob
   }
   N <- sum( obj$inputDataFreqInt )
   if( type == "prob" ) {
      fit <- prob
   } else {
      N <- sum( obj$inputDataFreqInt )
      mu <- N * prob
      if( type == "mean" ) {
         fit <- mu
      } else {
         eta <- numeric( length( NROW(obj$mfTrue) ) )
         eta[ obj$strZero ] <- as.double(-Inf)
         eta[ !obj$strZero ] <- log( mu[ !obj$strZero ] )
         fit <- eta
      }
   }
   if( mfTrue ) {
      result <- obj$mfTrue
      result$fit <- fit
   } else result <- fit
   result
}

get.imputedFreq <- function( obj, msgIfNone = TRUE ){
   stopifnot( inherits(obj, "cvam") )
   if( obj$method == "EM" ) {
      if( msgIfNone ) message( gettext(
         "No imputations are produced when method = 'EM'" ),
         domain = NA )
      return( invisible() )
   }
   if( obj$method == "MCMC" && obj$control$imputeEvery == 0L ) {
      if( msgIfNone ) message( gettext(
         "No imputations are produced when control$imputeEvery = 0" ),
         domain = NA )
      return( invisible() )
   }
   if( obj$method == "approxBayes" && ! obj$control$imputeApproxBayes ) {
      if( msgIfNone ) message( gettext(
      "No imputations are produced when control$imputeApproxBayes is FALSE" ),
         domain = NA )
      return( invisible() )
   }
   if( obj$nImpActual < 1L ) {
      if( msgIfNone ) message( gettext(
         "Insufficient iterations after burn-in, no imputations available" ),
         domain = NA )
      return( invisible() )
   }
   result <- obj$mfTrue
   result$freq <- NULL
   imp <- obj$imputedFreqInt
   colnames(imp) <- paste( "imp", as.character(1:NCOL(imp)), sep=".")
   result <- cbind( result, imp )
   result
}

get.minus2logPSeries <- function( obj, startValShift = TRUE,
   msgIfNone = TRUE, coda = ( obj$method == "MCMC" ) ){
   stopifnot( inherits(obj, "cvam") )
   if( obj$method == "EM" ) {
      if( msgIfNone ) message( gettext(
         "No series are produced when method = 'EM'" ),
         domain = NA )
      return( invisible() )
   }
   if( obj$nSampleActual < 1L ) {
      if( msgIfNone ) message( gettext(
         "Insufficient itertions after burn-in, no series is available" ),
         domain = NA )
      return( invisible() )
   }
   shift <- if( startValShift && obj$beganAtMode ) obj$startLogP else 0
   result <- -2 * ( obj$logPSeries - shift )
   if(coda) coda::mcmc( result, thin=obj$control$thinMCMC ) else result
}


get.coefSeries <- function( obj, msgIfNone = TRUE,
   coda = ( obj$method == "MCMC" ) ){
   stopifnot( inherits(obj, "cvam") )
   if( attr(obj$model, "saturated") ) {
      if( msgIfNone ) message( gettext(
         "Model is saturated, coefficients are not defined" ), domain = NA )
      return( invisible() )
   }
   if( obj$method == "EM" ) {
      if( msgIfNone ) message( gettext(
         "No series are stored when method = 'EM'" ),
         domain = NA )
      return( invisible() )
   }
   if( obj$nSampleActual < 1L ) {
      if( msgIfNone ) message( gettext(
         "Insufficient iterations after burn-in, no series available" ),
         domain = NA )
      return( invisible() )
   }
   result <- obj$betaSeries
   if(coda) coda::mcmc( result, thin=obj$control$thinMCMC ) else result
}


get.probSeries <- function( obj, levelNames=TRUE, sep=".",
   msgIfNone = TRUE, coda = ( obj$method == "MCMC" ) ){
   stopifnot( inherits(obj, "cvam") )
   if( obj$method == "EM" ) {
      if( msgIfNone ) message( gettext(
         "No series are stored when method = 'EM'" ),
         domain = NA )
      return( invisible() )
   }
   if( ! obj$control$saveProbSeries ) {
      if( msgIfNone ) message( gettext(
         "No probSeries was stored; control$saveProbSeries is FALSE" ),
         domain = NA )
      return( invisible() )
   }
   if( obj$nSampleActual < 1L ) {
      if( msgIfNone ) message( gettext(
         "Insufficient iterations after burn-in, no series available" ),
         domain = NA )
      return( invisible() )
   }
   result <- obj$probSeries
   if( levelNames ){
      sep <- as.character(sep)[1L]
      mfTrue <- obj$mfTrue
      mfTrue$freq <- NULL
      colnames(result)  <- interaction(mfTrue, sep=sep)
   }
   if(coda) coda::mcmc( result, thin=obj$control$thinMCMC ) else result
}

