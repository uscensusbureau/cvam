###########################################################################
# Unit tests using the UCBAdmissions data
###########################################################################
library(cvam)
dF <- microUCBAdmissions  # to save typing
# this reproduces the 3-way table UCBAdmissions
result <- table( Admit = dF$Admit, 
   Gender = dF$Gender, Dept = dF$Dept )
str(result)
print( all.equal( result, UCBAdmissions ) )

# do the same thing with xtabs, which accepts formula notation
result <- xtabs( ~ Admit + Gender + Dept, data=microUCBAdmissions )
print( all( result==UCBAdmissions ) )

# create a Freq variable and fill it with ones
microUCBAdmissions$Freq <- 1
# use aggregate to sum the Freq variable within categories of
# Admit, Gender, and Dept
result <- aggregate( Freq ~ Admit + Gender + Dept,
   data=microUCBAdmissions, FUN=sum )
print( head(result) )


# fit with glm
dF <- as.data.frame(UCBAdmissions)
M0 <- glm( Freq ~ Dept*Gender + Dept*Admit, family=poisson(), data=dF )
M1 <- glm( Freq ~ Dept*Gender + Dept*Admit + Gender*Admit,
   family=poisson(), data=dF )
M2 <- glm( Freq ~ Dept*Gender*Admit, family=poisson(), data=dF )
dF$muHat0 <- predict(M0, type="response")
dF$muHat1 <- predict(M1, type="response")
dF$muHat2 <- predict(M2, type="response")
fit0 <- xtabs( muHat0 ~ Admit + Gender + Dept, data=dF )
fit1 <- xtabs( muHat1 ~ Admit + Gender + Dept, data=dF )
fit2 <- xtabs( muHat2 ~ Admit + Gender + Dept, data=dF )
# under M0, the fitted conditional OR's should be 1.0:
print( fit0[1,1,] * fit0[2,2,] / ( fit0[1,2,] * fit0[2,1,] ) )
# under M1, the fitted conditional OR's should be equal:
print( fit1[1,1,] * fit1[2,2,] / ( fit1[1,2,] * fit1[2,1,] ) )
# under M2, the fitted conditional OR's should vary, and they
# should agree with corresponding OR's based on the observed
# frequencies, because M2 is saturated:
print( fit2[1,1,] * fit2[2,2,] / ( fit2[1,2,] * fit2[2,1,] ) )
print( anova(M0,M1,M2) )
d01 <- deviance(M0)-deviance(M1)
d12 <- deviance(M1)-deviance(M2)
d02 <- deviance(M0)-deviance(M2)
# make a list of 6 data frames, one per department
list2x2 <- as.list(1:6)
for( j in 1:6 ) list2x2[[j]] <- subset(dF, Dept==levels(dF$Dept)[j]  )
# function for computing deviance for LR test of independence
# within a department
myFunc <- function( dF ) {
   M <- glm( Freq ~ Gender + Admit, family=poisson(), data=dF )
   deviance(M)
}
# apply LR test to each department, returning a vector of deviances
dev <- sapply( list2x2, myFunc )
print( dev )
print( sum(dev) )
# should be same
print( anova(M0,M2) )

# run same models in cvam
dF <- as.data.frame(UCBAdmissions)
M0 <- cvam( ~ Dept*Gender + Dept*Admit, data=dF, freq=Freq )
M1 <- cvam( ~ Dept*Gender + Dept*Admit + Gender*Admit, data=dF, freq=Freq )
M2 <- cvam( ~ Dept*Gender*Admit, data=dF, freq=Freq )
anova(M0,M1,M2)
print( get.coef(M0, withSE=TRUE) )
print( head( get.fitted(M0, type="mean" ) ) )

# refit M0 with microdata to see that results are the same
M0 <- cvam( ~ Dept*Gender + Dept*Admit, data=microUCBAdmissions )
print( get.coef(M0, withSE=TRUE) )

# compute the deviance for model M0
M0 <- cvam( ~ Dept*Gender + Dept*Admit, data=dF, freq=Freq )
M2 <- cvam( ~ Dept*Gender*Admit, data=dF, freq=Freq, saturated=TRUE )
dev.M0 <- -2 * ( get.loglik(M0) - get.loglik(M2) )
print( dev.M0 )

# fit M1 as a conditional model
M1 <- cvam( ~ Dept*Gender + Dept*Admit + Gender*Admit | Dept + Gender,
   data=dF, freq=Freq )
print( head( get.fitted(M1, type="prob") ) )
