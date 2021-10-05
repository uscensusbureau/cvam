###########################################################################
# Unit tests using the abortion2000 data
###########################################################################
library(cvam)
# Create the coarsened factor RH, as described in the vignette
CenRace <- abortion2000$CenRace
Hisp <- abortion2000$Hisp
CenRace <- addNA(CenRace)
Hisp <- addNA(Hisp)
RH <- Hisp:CenRace
RH <- droplevels(RH)
levels(RH) <- list(
   nonHispWhite = "nonHisp:White",
   nonHispBlack = "nonHisp:Black",
   nonHispOther = "nonHisp:Other",
   Hisp = c("Hisp:White", "Hisp:Black", "Hisp:Hisp", "Hisp:Other", "Hisp:NA"),
   nonHispNA = "nonHisp:NA",
   NAWhite = "NA:White" )
RH  <- coarsened( RH, levelsList = list(
   nonHispNA = c("nonHispWhite", "nonHispBlack", "nonHispOther"), 
   NAWhite = c("nonHispWhite", "Hisp" ) ) )
abortion2000$RH <- RH

# plain table
print( summary(RH) )

# functions for extracting attributes
print( baseLevels(RH) )
print( nBaseLevels(RH) )
print( coarseLevels(RH) )
print( nCoarseLevels(RH) )
print( mapping(RH) )
print( is.latentFactor(RH) )

# is.naCoarsened
print( table(is.naCoarsened(RH)) )
RH2 <- RH
RH2[1:10] <- NA
print( table(is.naCoarsened(RH2)) )


# put in a data frame, make sure that attributes are preserved
dF <- data.frame( RH, freq=1 )
print( head(dF) )
print( identical( attributes(RH), attributes(dF$RH) ) ) 

# table, xtabs, ftable
PolViews <- abortion2000$PolViews
# there are some missing values
print( table( is.na(PolViews) ) )
# but the NAs don't show up in a table
print( table(RH, PolViews) )
# show the NAs
print( table(RH, PolViews, exclude=NULL) )
print( table(RH, PolViews = addNA(PolViews) ) )
print( table(RH, PolViews = coarsened(PolViews) ) )
# with xtabs
print( xtabs( ~ RH + PolViews, addNA=TRUE) )
# ftable
Sex <- abortion2000$Sex
print( ftable( addNA(PolViews) ~ Sex + RH, exclude=NULL ) )
# dropping the coarse levels
print( table( RH=dropCoarseLevels(RH), PolViews) )

# grouping data
groupedData = as.data.frame( xtabs( ~ Age + Sex + CenRace + Hisp,
   data=abortion2000, addNA=TRUE) )
groupedData <- subset( groupedData, Freq > 0 ) 
print( dim(groupedData) )
print( head(groupedData) )

CenRace <- addNA(groupedData$CenRace)
Hisp <- addNA(groupedData$Hisp)
RH <- Hisp:CenRace
RH <- droplevels(RH)
levels(RH) <- list(
   nonHispWhite = "nonHisp:White",
   nonHispBlack = "nonHisp:Black",
   nonHispOther = "nonHisp:Other",
   Hisp = c("Hisp:White", "Hisp:Black", "Hisp:Hisp", "Hisp:Other", "Hisp:NA"),
   nonHispNA = "nonHisp:NA",
   NAWhite = "NA:White" )
RH  <- coarsened( RH, levelsList = list(
   nonHispNA = c("nonHispWhite", "nonHispBlack", "nonHispOther"), 
   NAWhite = c("nonHispWhite", "Hisp" ) ) )
groupedData$RH <- RH
print( aggregate( Freq ~ RH, FUN=sum, data=groupedData) )

# aggregate
dFaggr <- aggregate( freq ~ RH, FUN=sum, data=dF )
print( identical( attributes(dF$RH), attributes(dFaggr$RH) ) ) 

# quickly written function that accepts a single coarsened factor,
# optionally with frequencies, and computes ML estimates for the
# base-level probabilities. This function is used for unit testing.
quickEM <- function( obj, freq = rep(1, length(obj)),
   maxits=1000, eps=1e-06 ){
   # identify the name of the coarsened factor
   mc <- match.call()
   objName <- as.character( mc[[2]] )
   if( objName == "freq" ) stop( gettext(
      "Main argument cannot be named 'freq'"), domain = NA )
   # check args
   stopifnot( inherits(obj, "coarsened") )
   stopifnot( length(freq)==length(obj) )
   stopifnot( all( !is.na(freq) ) )
   stopifnot( all(freq>=0) )
   stopifnot( maxits > 0 )
   # aggregate the data to create model frame, then 
   # pull out the coarsened factor and frequencies
#   obj <- sticky::sticky(obj)
   mf <- aggregate( freq ~ obj, FUN=sum )
   names(mf)[ names(mf)=="obj" ] <- objName
   cFac <- mf[[objName]]
   freq <- mf[["freq"]]
   # starting value: uniform probabilities
   theta <- rep( 1 / nBaseLevels(cFac), nBaseLevels(cFac) )
   names(theta) <- baseLevels(cFac)
   # prepare for iteration
   Ex <- theta.new <- theta
   iter <- 0
   converged <- FALSE
   llvec <- numeric(maxits)
   llvec[] <- NA
   while( ( ! converged ) & ( iter < maxits ) ) {
      iter <- iter + 1
      theta <- theta.new
      Ex[] <- loglik <- 0
      # e-step
      for( i in seq_along(cFac) ){
         iC <- as.integer( cFac[i] )
         if( iC %in% attr(cFac,"baseLevelCodes") ) {
            Ex[iC] <- Ex[iC] + freq[i]
            loglik <- loglik + freq[i] * log( theta[iC] )
         } else {
            w <- match(iC, attr(cFac, "coarseLevelCodes") )
            w <- ( mapping(cFac)[w,] == 1 )
            Ex[w] <- Ex[w] + freq[i] * theta[w] / sum( theta[w] )
            loglik <- loglik + freq[i] * log( sum(theta[w]) )
         }
      }
      llvec[iter] <- loglik
      names(loglik) <- NULL
      # m-step
      theta.new <- Ex / sum(Ex)
      # convergence check
      converged <- all( abs(theta.new-theta) <= eps*abs(theta) )      
   }
   list( theta=theta, iter=iter, converged=converged, 
      llvec=llvec[1:iter], loglik=loglik )
}

# run quickEM on microdata
RH <- abortion2000$RH
resultA <- quickEM(RH, eps=1e-08)
print(resultA)

# run cvam on microdata
resultB <- cvam( ~ RH, est = ~ RH, data=abortion2000 )
print( summary(resultB, showEst=TRUE) )

# run quickEM on grouped data
RH <- groupedData$RH
Freq <- groupedData$Freq
resultC <- quickEM(RH, Freq, eps=1e-08)
print(resultC)

# run cvam on grouped data
resultD <- cvam( ~ RH, est = ~ RH, data=groupedData, freq=Freq )
print( summary(resultD) )



## four-way model
CenRace <- addNA(abortion2000$CenRace)
Hisp <- addNA(abortion2000$Hisp)
RH <- Hisp:CenRace
RH <- droplevels(RH)
levels(RH) <- list(
   nonHispWhite = "nonHisp:White",
   nonHispBlack = "nonHisp:Black",
   nonHispOther = "nonHisp:Other",
   Hisp = c("Hisp:White", "Hisp:Black", "Hisp:Hisp", "Hisp:Other", "Hisp:NA"),
   nonHispNA = "nonHisp:NA",
   NAWhite = "NA:White" )
RH  <- coarsened( RH, levelsList = list(
   nonHispNA = c("nonHispWhite", "nonHispBlack", "nonHispOther"), 
   NAWhite = c("nonHispWhite", "Hisp" ) ) )

# copy the four variables into a data frame
dF <- data.frame( Sex = abortion2000$Sex, RH = RH,
   PolViews = abortion2000$PolViews, AbAny = abortion2000$AbAny )

myFormula <- ~ Sex*RH*PolViews + AbAny*Sex + AbAny*RH + AbAny*PolViews
myMod <- cvam( myFormula, data=dF )
satMod <- cvam( ~ Sex*RH*PolViews*AbAny, data=dF, saturated=TRUE )
print( anova( myMod, satMod, pval=TRUE ) )
print( anova( myMod, satMod, method="AIC" ) )
# compute and summarize the fitted values
muHat <- get.fitted(myMod, type="mean")$fit
print( summary( muHat ) )

noSex <- cvam( ~ Sex*RH*PolViews + AbAny*RH + AbAny*PolViews, data=dF)
print( anova( noSex, myMod, pval=TRUE ) )
noRH  <- cvam( ~ Sex*RH*PolViews + AbAny*Sex + AbAny*PolViews, data=dF)
print( anova( noRH, myMod, pval=TRUE ) )
noPol <- cvam( ~ Sex*RH*PolViews + AbAny*Sex + AbAny*RH, data=dF)
print( anova( noPol, myMod, pval=TRUE ) )

# use a boundary criterion that is less strict
satMod <- cvam( ~ Sex*RH*PolViews*AbAny, data=dF, saturated=TRUE,
   control=list(critBoundary=1e+06 ) )
print( round( get.fitted(satMod, type="prob")$fit, 6) )

myPrior <- cvamPrior( flatten=7.2 )
# re-fit and compare models using the flattening constant
myMod <- cvam( myFormula, data=dF, prior=myPrior )
satMod <- cvam( ~ Sex*RH*PolViews*AbAny, data=dF,
   saturated=TRUE, prior=myPrior )
print( anova( myMod, satMod, pval=TRUE, method="logP") )

# re-run EM with a ridge factor of 0.5
fitEM <- cvam( ~ Sex * RH * PolViews * AbAny, data=dF )
fitEM.ridge <- cvam( ~ Sex * RH * PolViews * AbAny, data=dF,
   prior=cvamPrior( ridge=.5 ) )
round( get.fitted(fitEM.ridge, type="prob", mfTrue=FALSE ), 5)
print( head( get.coef(fitEM.ridge, withSE=TRUE) ) )
print( -2 * ( get.loglik(fitEM.ridge) - get.loglik(fitEM) ) )
print( exp( get.loglik(fitEM) - get.loglik(fitEM.ridge) ) )

set.seed(87900)
fitMCMC <- cvam( fitEM.ridge, method="MCMC" )

set.seed(87900)
fitMCMC <- cvam( fitEM.ridge, method="MCMC",
   control=list( typeMCMC="RWM", tuneRWM=c(1000,.17) ) )
print( summary(fitMCMC) )
