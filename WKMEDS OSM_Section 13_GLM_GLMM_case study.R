#############################################################################
##### ICES Workshop on Methods for Estimating Discard Survival (WKMEDS) #####
#####  Online Supporting Materials for the Cooperative Research Report  #####
#####                                                                   #####
#####       Section 13.1 - Case study to illustrate the use of          #####
#####                 binomial GLM and GLMM models                      #####
#####                                                                   #####
#####                 Version 2.0      1st July 2021                    #####
#############################################################################

# Original code by H Benoît 2018-08-24
# Edited by M Breen 2021-07-01

# This code supports the examples and illustrations presented in section 13.1 of 
# the WKMEDS Cooperative Research Report (CRR).


#### Data SOurce
# Data are for American plaice in moribund condition prior to release in water
# from the study of Benoît et al. 2012. Estimating fishery-scale rates of
# discard mortality using conditional reasoning. Fish. Res. 125-126: 318-330.

### -------------------------------------------------------------------- ###

## Prep & Required Packages

rm(list=ls())

library(lme4)


## Set Working Directory
# This code is written assumming the user has it and  the supporting data in a working directory.

setwd("E:/IMR/UM/ICES WKMEDS/CRR Report/WKMEDS CRR Final Edit/Online Support Materials/Chapter_13_Modelling Survival Data/GLM and GLMM examples")

### -------------------------------------------------------------------- ###


#!=========================================================================!#
#!=======                     Section 13.1.1                         ======!#
#!=======  13.1.1	Generalized linear model (GLM) for survival data   ======!#
#!=========================================================================!#

#!=========================================================================!#
#!======= Info Box 13.1. Fitting a GLM in R - worked example          =====!#
#!=========================================================================!#


### Read in data
pla2=read.csv("Benoit et al_American plaice.csv")
View(pla2)

### Fit the basic binomial GLM to the American plaice data
glm1=glm(dead~decktime,family=binomial(link='logit') ,data=pla2)
summary(glm1)

## look at residuals and consider evidence for over-dispersion
E1 <- resid(glm1, type = "pearson")
plot(x = pla2$decktime, y = E1)
N  <- nrow(pla2)
p  <- length(coef(glm1))
sum(E1^2) / (N - p)


## Generate predicted values and confidence intervals
df = data.frame( decktime=rep((dt=0:75),times=1)) 
preds = predict(glm1, newdata = df, type = "link", se.fit = TRUE)
critval = 1.96 ## approx 95% CI
upr = preds$fit + (critval * preds$se.fit)
lwr = preds$fit - (critval * preds$se.fit)
fit = preds$fit
fit2 = glm1$family$linkinv(fit)
upr2 = glm1$family$linkinv(upr)
lwr2 = glm1$family$linkinv(lwr)

## Summarize the observations
pla2$dt=10*round(pla2$decktime/10)
num=as.data.frame(aggregate(list(num=pla2$dead),by=list(dt=pla2$dt),sum))
den=as.data.frame(aggregate(list(den=pla2$dead),by=list(dt=pla2$dt),length))
prop=merge(num,den)
prop$prop=prop$num/prop$den

### Plot the predictions and summarized observations (Figure 13.1)

## Open separate graphics window
# windows()

## Set up Tiff file export
# to working directory
# tiff("binomial GLM_decktime.tif", width = 4, height = 4, units = "in", res = 600, pointsize = 10, compression="lzw")

with(prop,plot(dt,prop,cex=log10(den),xlim=c(0,75),ylim=c(0,1),xaxs="i",yaxs="i",xlab='Deck time (min)',ylab='Proportion dead at time of discarding'))
polygon(c(0:75, 75:0), c(lwr2[1:76], upr2[76:1]), col =  rgb(0, 0, 0, alpha = 0.1), border = NA)
polygon(c(0:75, 75:0), c(lwr2[77:152], upr2[152:77]), col =  rgb(1, 0, 0, alpha = 0.1), border = NA)
lines(dt,fit2[1:76])
lines(dt,fit2[77:152],col='red')
legend(50,0.3,c('1','10','100'),pch=rep(1,3),pt.cex=log10(c(1.1,10,100)),bty='n',cex=0.8)

##Close Tiff file export
# dev.off()

### -------------------------------------------------------------------- ###


#!=========================================================================!#
#!======= Info Box 13.2. Fitting a GLM in R - worked example part 2   =====!#
#!=========================================================================!#


### Consider other other models, using forward selection based on AIC
glm2=glm(dead~length,family=binomial(link='logit') ,data=pla2)

glm3=glm(dead~season,family=binomial(link='logit') ,data=pla2)

AIC(glm1,glm2,glm3) #the model with decktime provides the best of the single covariate models based on AIC 

glm4=glm(dead~decktime+length,family=binomial(link='logit') ,data=pla2)

glm5=glm(dead~decktime+season,family=binomial(link='logit') ,data=pla2)

AIC(glm1,glm4,glm5) #there is much more evidence for a model with decktime and season compared to glm1 
summary(glm5)

glm6=glm(dead~decktime+season+length,family=binomial(link='logit') ,data=pla2)
summary(glm6)

AIC(glm5,glm6) #the addition of body length does not significantly improve the model
anova(glm5,glm6,test='Chi')

AIC(glm1,glm2,glm3,glm4,glm5,glm6)


## look at residuals and consider evidence for over-dispersion
E5 <- resid(glm5, type = "pearson")
N  <- nrow(pla2)
p  <- length(coef(glm5))
sum(E5^2) / (N - p)

## plot the predictions and summarized observations

pla2$colour='black'
pla2$colour[pla2$season=='summer']='darkgrey'

# tiff("Figure 13_4.tif", width = 4, height = 4, units = "in", res = 600, pointsize = 10, compression="lzw")

plot(x = pla2$decktime, y = E5,col=pla2$colour,xlab='Deck time (min)',ylab='Pearson residual (grey=summer)')

# dev.off()


#summarize observations
num=as.data.frame(aggregate(list(num=pla2$dead),by=list(dt=pla2$dt,season=pla2$season),sum))
den=as.data.frame(aggregate(list(den=pla2$dead),by=list(dt=pla2$dt,season=pla2$season),length))
prop=merge(num,den)
prop$prop=prop$num/prop$den

df = data.frame( decktime=rep((dt=0:75),times=2), season=c(rep('summer',76),rep('autumn',76))) 
preds = predict(glm5, newdata = df, type = "link", se.fit = TRUE)
critval = 1.96 ## approx 95% CI
upr = preds$fit + (critval * preds$se.fit)
lwr = preds$fit - (critval * preds$se.fit)
fit = preds$fit
fit2 = glm5$family$linkinv(fit)
upr2 = glm5$family$linkinv(upr)
lwr2 = glm5$family$linkinv(lwr)

# windows()
# tiff("Figure 13_3.tif", width = 4, height = 4, units = "in", res = 600, pointsize = 10, compression="lzw")

with(prop[prop$season=='summer',],plot(dt,prop,cex=log10(den),xlim=c(0,75),ylim=c(0,1),xaxs="i",yaxs="i",xlab='Deck time (min)',ylab='Proportion dead at time of discarding'))
with(prop[prop$season=='autumn',],points(dt,prop,cex=log10(den),col='darkgrey'))
polygon(c(0:75, 75:0), c(lwr2[1:76], upr2[76:1]), col =  rgb(0, 0, 0, alpha = 0.1), border = NA)
polygon(c(0:75, 75:0), c(lwr2[77:152], upr2[152:77]), col =  rgb(0.5, 0.5, 0.5, alpha = 0.1), border = NA)
lines(dt,fit2[1:76])
lines(dt,fit2[77:152],col='darkgrey')
legend(38,0.3,c('1','10','100'),pch=rep(1,3),pt.cex=log10(c(1.1,10,100)),bty='n',cex=0.8)
legend(50,0.3,c('Summer','Autumn'),lty=c(1,1),col=c('darkgrey','black'),bty='n',cex=0.8)

# dev.off()

### -------------------------------------------------------------------- ###


#!=========================================================================!#
#!=======                     Section 13.1.2                         ======!#
#!=======       13.1.2	Generalized linear mixed models (GLMM)       ======!#
#!=========================================================================!#

#!=========================================================================!#
#!======= Info Box 13.4. Fitting a GLMM in R - worked example part 3 ======!#
#!=========================================================================!#

### Required package
library(lme4)

### GLMM with random intercept
glmm1=glmer(dead~decktime+season + (1|haul) ,data=pla2,family=binomial )
summary(glmm1)

### GLMM with random slope
glmm2=glmer(dead~decktime+season + (-1+decktime|haul) ,data=pla2,family=binomial )
summary(glmm2)
AIC(glmm1,glmm2)

### GLMM with random slope & intercept
glmm3=glmer(dead~decktime+season + (decktime|haul),data=pla2,family=binomial )
summary(glmm3)
anova(glmm1,glmm2,glmm3)

### GLMM with random slope & intercept, without the effect of season
glmm4=glmer(dead~decktime + (decktime|haul),data=pla2,family=binomial )
summary(glmm4)
anova(glmm3,glmm4)

## plot the residuals
EM4 <- resid(glmm4, type = "pearson")

# tiff("binomial GLMM_resid.tif", width = 4, height = 4, units = "in", res = 600, pointsize = 10, compression="lzw")

plot(x = pla2$decktime, y = EM4,xlab='Deck time (min)',ylab='Pearson residual')

# dev.off()


### plot the predictions and summarized observations

## Create data on grid, and the matching X matrix
MyData =  data.frame( decktime=rep((dt=0:75),times=1)) 
X = model.matrix(~decktime , data = MyData)

## Extract parameters and parameter covariance matrix
FEs    = fixef(glmm4)
CovFEs = vcov(glmm4)

## Calculate the fitted values in the predictor scale
MyData$eta = X %*% FEs
MyData$Pi  = exp(MyData$eta) / (1 + exp(MyData$eta))

## Calculate the CIs on the scale of the predictor function
MyData$se   = sqrt(diag(X %*% CovFEs %*% t(X)))
MyData$CIUp  = exp(MyData$eta + 1.96 *MyData$se) / (1 + exp(MyData$eta  + 1.96 *MyData$se))
MyData$CILo  = exp(MyData$eta - 1.96 *MyData$se) / (1 + exp(MyData$eta  - 1.96 *MyData$se))

## Extract the random effects for plotting
REs=ranef(glmm4)$haul

# tiff("binomial GLMM.tif", width = 4, height = 4, units = "in", res = 600, pointsize = 10, compression="lzw")
# windows()

with(prop,plot(dt,prop,cex=log10(den),xlim=c(0,75),ylim=c(0,1),xaxs="i",yaxs="i",xlab='Deck time (min)',ylab='Proportion dead at time of discarding'))
polygon(c(0:75, 75:0), c(MyData$CILo[1:76], MyData$CIUp[76:1]), col =  rgb(0, 0, 0, alpha = 0.1), border = NA)

for(ii in 1:nrow(REs)){
  REval=REs[ii,]
  pred=1/(1+exp(-(REval[1,1]+FEs[1] + (REval[1,2]+FEs[2])*MyData$decktime  )))
  lines(MyData$decktime,pred,lwd=0.1,lty=3)
}
lines(dt,MyData$Pi,lwd=2)
legend(50,0.3,c('1','10','100'),pch=rep(1,3),pt.cex=log10(c(1.1,10,100)),bty='n',cex=0.8)

# dev.off()


### -------------------------------------------------------------------- ###
### -------------------------------------------------------------------- ###

