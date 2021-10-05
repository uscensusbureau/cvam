library(cvam)
hivtest$L <- latentFactor( NROW(hivtest), 2 )
set.seed(125)
fit <- cvam( ~ L*A + L*B + L*C + L*D, data=hivtest, freq=COUNT,
   control = list( startValJitter=.1 ) )
print( summary(fit) )

print( cvamEstimate( list( ~L, ~A|L, ~B|L, ~C|L, ~D|L ), fit ) )

# perform the lack-of-fit test
fitSat <- cvam( ~ A*B*C*D, data=hivtest, freq=COUNT )
print( anova( fit, fitSat, pval=TRUE ) )

satFrame <- get.fitted( fitSat, type="mean" )
# get rid of the fitted values, because they are redundant
satFrame$fit <- NULL
LCFrame <-  get.fitted( fit, type="mean" )
muHatTable <- xtabs( fit ~ A + B + C + D, data=LCFrame )
muHatFrame <- as.data.frame( muHatTable, responseName = "muHat" )
muHat <- muHatFrame$muHat
quasiPearson <- ( satFrame$freq - muHat ) / sqrt( muHat )
satFrame$muHat <- round( muHat, 3 )
satFrame$quasiPearson <- round( quasiPearson, 2 )
print( satFrame )

set.seed(85657)
fitLAB <- cvam( ~ L*A + L*B + L*C + L*D + L*A*B, 
   data=hivtest, freq=COUNT,
   control = list(startValJitter=.1) )
print( anova(fit, fitLAB, fitSat, pval=TRUE) )
fitLAC <- cvam( ~ L*A + L*B + L*C + L*D + L*A*C, 
   data=hivtest, freq=COUNT,
   control = list(startValJitter=.1) )
print( anova(fit, fitLAC, fitSat, pval=TRUE) )
fitLAD <- cvam( ~ L*A + L*B + L*C + L*D + L*A*D, 
   data=hivtest, freq=COUNT,
   control = list(startValJitter=.1) )
print( anova(fit, fitLAD, fitSat, pval=TRUE) )
fitLBC <- cvam( ~ L*A + L*B + L*C + L*D + L*B*C, 
   data=hivtest, freq=COUNT,
   control = list(startValJitter=.1) )
print( anova(fit, fitLBC, fitSat, pval=TRUE) )
fitLBD <- cvam( ~ L*A + L*B + L*C + L*D + L*B*D, 
   data=hivtest, freq=COUNT,
   control = list(startValJitter=.1) )
print( anova(fit, fitLBD, fitSat, pval=TRUE) )
fitLCD <- cvam( ~ L*A + L*B + L*C + L*D + L*C*D, 
   data=hivtest, freq=COUNT,
   control = list(startValJitter=.1) )
print( anova(fit, fitLCD, fitSat, pval=TRUE) )
fitBoth <- cvam( ~ L*A + L*B + L*C + L*D + L*A*D + L*B*C, 
   data=hivtest, freq=COUNT,
   control = list(startValJitter=.1) )
print( anova(fitLAD, fitBoth) )
print( anova(fitLBC, fitBoth) )

# get predicted probabilities and display them with the dataset
pred <- cvamPredict( ~L, fitLBC, data=hivtest )
print( cbind( hivtest, round(pred, 3) ) )

predFrame <- hivtest[1:8,]
predFrame$COUNT <- NULL
predFrame[["A"]][] <- NA
predFrame[["B"]][] <- NA
predFrame[["C"]][] <- NA
predFrame[["D"]][] <- NA
predFrame[["A"]][1] <- "pos"; predFrame[["A"]][2] <- "neg"
predFrame[["B"]][3] <- "pos"; predFrame[["B"]][4] <- "neg"
predFrame[["C"]][5] <- "pos"; predFrame[["C"]][6] <- "neg"
predFrame[["D"]][7] <- "pos"; predFrame[["D"]][8] <- "neg"
predFrame[["A"]] <- coarsened( predFrame[["A"]] )
predFrame[["B"]] <- coarsened( predFrame[["B"]] )
predFrame[["C"]] <- coarsened( predFrame[["C"]] )
predFrame[["D"]] <- coarsened( predFrame[["D"]] )
pred <- cvamPredict( ~L, fitLBC, data=predFrame )
print( cbind( predFrame, round(pred, 3) ) )

pred <- cvamPredict( ~L, fit, data=predFrame )
print( cbind( predFrame, round(pred, 3) ) )


# re-fit the model with EM using a small ridge factor
set.seed(7666)
fitLBC <- cvam( ~ L*A + L*B + L*C + L*D + L*B*C, 
   data=hivtest, freq=COUNT, prior=cvamPrior( ridge=.1 ),
   control = list(startValJitter=.1) )
# do a long run of MCMC and save ten imputed datasets
fitMCMC <- cvam(fitLBC, method="MCMC",
   control=list( typeMCMC="RWM", tuneRWM=c(1000,.5),
      iterMCMC=25000, imputeEvery=2500 ) )
# check to see if any label switching has occurred
impData <- get.imputedFreq(fitMCMC)
print( head(impData) )

impData$freq <- impData[["imp.1"]] # first imputation
BCL <- xtabs( freq ~ B + C + L, data=impData )
print( BCL )

# use multiple imputations to examine the conditional
# BC odds ratios given L=1 and L=2
est.list <- SE.list <- as.list(1:10)
for( m in 1:10 ) {
   # get the imputed marginal table BxCxL
   impName <- paste( "imp", format(m), sep="." )
   impData$freq <- impData[[impName]]
   BCL <- xtabs( freq ~ B + C + L, data=impData )
   # add 1/2 to every cell to avoid problems
   BCL <- BCL + .5
   # get BC log-odds ratio and SE for L=1
   BCL.1 <- BCL[,,"1"]
   logOR.1 <- log( ( BCL.1[1,1] * BCL.1[2,2] ) /
      ( BCL.1[1,2] * BCL.1[2,1] ) )
   SE.1 <- sqrt( sum( 1/BCL.1 ) )
   # get BC log-odds ratio and SE for L=2
   BCL.2 <- BCL[,,"2"]
   logOR.2 <- log( ( BCL.2[1,1] * BCL.2[2,2] ) /
      ( BCL.2[1,2] * BCL.2[2,1] ) )
   SE.2 <- sqrt( sum( 1/BCL.2 ) )
   # save the estimates and SEs
   est.list[[m]] <- c( logOR.1, logOR.2 )
   SE.list[[m]] <- c( SE.1, SE.2 )
}
print( miInference( est.list, SE.list ) )
