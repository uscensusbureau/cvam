library(cvam)
# fit the model of independence
M0 <- cvam( ~ V1 + V2, freq=n, data=crime )
# fit the model of non-independence
M1 <- cvam( ~ V1 * V2, freq=n, data=crime )
# compare them
print( anova(M0,M1, pval=TRUE) )

print( summary(M0) )

dF <- get.fitted(M0, type="mean")
print( dF )
print( ( dF$fit[1] * dF$fit[4] ) / ( dF$fit[2] * dF$fit[3] ) )

# compute the quasi-Pearson residuals
muHat <- get.fitted(M0, type="mean")$fit
fHatSat <- get.fitted(M1, type="mean")$freq
quasiPearson <- ( fHatSat - muHat ) / sqrt( muHat )
print( quasiPearson )

# fit the saturated model to the crime data
result <- cvam( ~ V1 * V2, data=crime, freq=n)
# run it again, starting from the previous result
result <- cvam(result)
print( summary(result, showCoef=FALSE) )

# fit the model of non-independence
fit <- cvam( ~ V1 * V2, data=crime, freq=n )
# display predictions for V1
print( cvamPredict( ~ V1, fit, data=crime ) )
# display predicted frequencies for V1
print( cvamPredict( ~ V1, fit, data=crime, freq=n ) )

# display predicted frequencies for V1 and V2
print( cvamPredict( ~ V1 + V2, fit, data=crime, freq=n ) )

set.seed(69852)
print( cvamImpute( fit, data=crime ) )
print( cvamImpute( fit, data=crime, freq=n ) )

# fit the non-independence model to the crime data
fitML <- cvam( ~ V1 * V2, data=crime, freq=n )
# display the ML estimate for beta and pi
print( get.coef( fitML ) )
print( get.fitted( fitML, type="prob" )$fit )
# draw from the approximate posterior, display new beta and pi
set.seed(83425)
obj <- cvam(fitML, method="approxBayes")
print( get.coef( obj ) )
print( get.fitted( obj, type="prob" )$fit )

obj <- cvam(fitML, method="approxBayes",
   control=list(iterApproxBayes=5000, saveProbSeries=TRUE) )
# display the first few beta and pi vectors
print( head( get.coefSeries(obj) ) )
print( head( get.probSeries(obj) ) )

pi.series <- get.probSeries(obj)
delta <- pi.series[,3] - pi.series[,2]
print( summary(delta) )
print( sum( delta > 0 ) )

set.seed(4358)
fit <- cvam( ~ V1 * V2, data=crime, freq=n, method="MCMC")
print( summary(fit) )

betaSeries <- get.coefSeries( fit )
print( summary( betaSeries ) )

print( get.fitted(fit, type="prob") )

set.seed(4358)
fit <- cvam( ~ V1 * V2, data=crime, freq=n, method="MCMC",
   control=list( saveProbSeries=TRUE ) )
piSeries <- get.probSeries(fit)
delta <- piSeries[,3] - piSeries[,2]
print( summary(delta) )
print( sum( delta > 0 ) )


impList <- as.list(1:10) # a list to store the imputed datasets
set.seed(769090)         # for reproducibility
for(m in 1:10) {
   # run MCMC under the non-independence model
   tmp <- cvam( ~ V1 * V2, data=crime, freq=n, method="MCMC")
   # impute under the simulated parameters
   impList[[m]] <- cvamImpute( tmp, crime, freq=n)
}
# display the first two imputations
print( impList[1:2] )

result <- cvam( ~ V1 * V2, data=crime, freq=n, method="MCMC",
   control=list( iterMCMC=5000, imputeEvery=500 ) )
print( get.imputedFreq(result) )

#  run EM, then create ten imputations with approxBayes
fitML <- cvam( ~ V1 * V2, data=crime, freq=n ) 
result <- cvam( fitML, method="approxBayes",
   control=list( iterApproxBayes=10, imputeApproxBayes=TRUE ) )
print( get.imputedFreq(result) )

set.seed(54981)
result <- cvam( fitML, method="MCMC",
   control=list( iterMCMC=5000, imputeEvery=500 ) )
impData <- get.imputedFreq(result)[-(1:2)] # just the frequencies 
est.list <- std.err.list <- as.list(1:10)  # to hold the estimates and SEs
for( m in 1:10 ) {
   f <- impData[,m]
   est.list[[m]] <- log( (f[1] * f[4]) / (f[2] * f[3]) )
   std.err.list[[m]] <- sqrt( sum(1/f) )
}
print( miInference( est.list, std.err.list ) )
