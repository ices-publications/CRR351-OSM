#############################################################################
##### ICES Workshop on Methods for Estimating Discard Survival (WKMEDS) #####
#####  Online Supporting Materials for the Cooperative Research Report  #####
#####                                                                   #####
#####           Section 13.2 - Analysis of Longitudinal Data            #####
#####                                                                   #####
#####                 Version 2.1       1st July 2021                   #####
#############################################################################

## Original code by H Benoît 2014-12-08, revised 2016-06-07, 2020-12-17
## Edited by M Breen 2021-06-24 and 2012-07-01

# This code supports the examples presented in section 13.2 of 
# the WKMEDS Cooperative Research Report (CRR).

### -------------------------------------------------------------------- ###


## Prep & Required Packages

rm(list=ls())
library(survival)
library(MASS)

## Set Working Directory
# This code is written assumming the user is running it and the supporting data from a working directory.

setwd("E:/IMR/UM/ICES WKMEDS/CRR Report/WKMEDS CRR Final Edit/Online Support Materials/Chapter_13_Modelling Survival Data/Longitudinal data examples")


#!=========================================================================!#
#!=======                     Section 13.2.2                         ======!#
#!=======     Non-parametric method: Kaplan-Meier survival model     ======!#
#!=========================================================================!#

#!=========================================================================!#
#!====    Info Box 13.5 - Kaplan-Meier (KM)  model - Worked example   =====!#
#!=========================================================================!#

# This worked example presents Kaplan-Meier analysis of the mortality data for 
# cod in excellent and poor condition (i.e. vitality), from Benoît et al. (2012) 

### Read in data
cod=read.table("sGSL_cod.csv", header=TRUE,sep=",")

## select the results for the two condition codes of interest
codKM=cod[cod$condition %in% c(1,3),] 

### Estimate the Kaplan-Meier survival functions 
km_est <- survfit(Surv(hours, (1-censored)) ~ condition, data=codKM, type="kaplan-meier") 
# note - the function uses an indicator for survival rather than for censoring 
# hence the use of (1-censored)


## log-rank test for differences between groups
survdiff(Surv(hours, (1-censored)) ~ condition, data=codKM)

## Wilcoxon test for differences between groups
survdiff(Surv(hours, (1-censored)) ~ condition, data=codKM,rho=1)


### -------------------------------------------------------------------- ###

### Plot the KM survival functions

# windows()
# tiff("Figure 13_7.tif", width = 4, height = 4, units = "in", res = 600, pointsize = 10, compression="lzw")

plot(km_est, main="Kaplan-Meier estimate for Excellent (black) and Poor (grey) vitality classes",lty=2,cex=1,cex.main=0.5,
     xlab="Hours", ylab="Survival function", col= c("black","darkgrey"))

## Plot the confidence intervals
# first create a matrix to index the correct rows to use from the km_est object 
rm=matrix(nrow=2,ncol=2)
rm[1,1]=1
rm[1,2]=as.numeric(km_est$strata[1])
rm[2,1]=rm[1,2]+1
rm[2,2]=rm[1,2]+as.numeric(km_est$strata[2])

censor=km_est$n.censor
censor[censor>1]=1

## get the data from those rows and plot as a semi-transparent polygon
cond=1
km_est1=as.data.frame(cbind(time=km_est$time[rm[cond,1]:rm[cond,2]],upper=km_est$upper[rm[cond,1]:rm[cond,2]],lower=km_est$lower[rm[cond,1]:rm[cond,2]]))
with(km_est1,polygon(c(time, time[length(time):1]), c(lower, upper[length(time):1]), col =  rgb(0, 0, 0, alpha = 0.35), border = NA))
points(km_est$time[c(rm[cond,1]:rm[cond,2])],km_est$surv[c(rm[cond,1]:rm[cond,2])],cex=censor[c(rm[cond,1]:rm[cond,2])],pch=3,col='black')

cond=2
km_est2=as.data.frame(cbind(time=km_est$time[rm[cond,1]:rm[cond,2]],upper=km_est$upper[rm[cond,1]:rm[cond,2]],lower=km_est$lower[rm[cond,1]:rm[cond,2]]))
with(km_est2,polygon(c(time, time[length(time):1]), c(lower, upper[length(time):1]), col =  rgb(0.5, 0.5, 0.5, alpha = 0.35), border = NA))
points(km_est$time[c(rm[cond,1]:rm[cond,2])],km_est$surv[c(rm[cond,1]:rm[cond,2])],cex=censor[c(rm[cond,1]:rm[cond,2])],pch=3,col='darkgrey')

# dev.off()

### -------------------------------------------------------------------- ###


#!=========================================================================!#
#!=======   Semi-parametric method: Cox proportional hazards model   ======!#
#!=======                     Section 13.2.3                         ======!#
#!=========================================================================!#

#!=========================================================================!#
#!====    Info Box 13.6 - The Cox proportional hazards (CPH) model    =====!#
#!====                         worked example                         =====!#
#!=========================================================================!#


# Cox proportional-hazards analysis of the mortality data for cod 
# in excellent and poor condition (i.e. vitality) from Benoît et al. (2012) 

### Estimate the predicted survival functions for the two condition classes
coxfit=coxph(Surv(hours, (1-censored)) ~ condition, data=codKM) 
summary(coxfit)

## test proportional hazards assumption
cox.zph(coxfit)

# => not significant p=0.057, therefore the proportional hazards assumption is valid.

### -------------------------------------------------------------------- ###

### Plot the CPH survival functions

# windows()
# tiff("Figure 13_8.tif", width = 4, height = 4, units = "in", res = 600, pointsize = 10, compression="lzw")

coxpred1=survfit(coxfit, newdata=data.frame(condition=3))
plot(coxpred1,col="darkgrey",main="Cox prop. hazards model estimate for Excellent (black) and Poor (grey) vitality classes",cex.main=0.5,xlab="Hours", ylab="Survival function")
censor=coxpred1$n.censor
censor[censor>1]=1
points(coxpred1$time,coxpred1$surv,cex=censor,pch=3,col='darkgrey')

coxpred2=survfit(coxfit, newdata=data.frame(condition=1))
lines(coxpred2,col="black")
censor=coxpred2$n.censor
censor[censor>1]=1
points(coxpred2$time,coxpred2$surv,cex=censor,pch=3,col='black')

# dev.off()

### -------------------------------------------------------------------- ###


#!=========================================================================!#
#!=======               Parametric modelling methods                 ======!#
#!=======                     Section 13.2.4                         ======!#
#!=========================================================================!#

#!=========================================================================!#
#!====       Info Box 13.7 - The Weibull model - Worked example       =====!#
#!=========================================================================!#


### estimate the model parameters
(sr=survreg(Surv(hours, (1-censored)) ~ condition, data=codKM,dist="weibull"))
sr$coefficients
sr$scale

## calculate AICc
AIC=-sr$loglik[2]*2 + 2*sr$df
AICc=AIC + (2*sr$df*(sr$df+1))/(sr$df.residual-1)


### use the estimates to plot the survival functions

## Weibull cumulative density function 
dweibull=function(x, theta, omega){
  G=1 - exp(-((theta*x)^omega))
  return(G)
}

### Plot the results
xx=seq(0, 110, len = 110)
theta1=exp(-sr$coefficients[1])
theta2=exp(-(sr$coefficients[1]+sr$coefficients[2]))
omega=sr$scale

## Plot curves for each injury code. 
# NOTE - this model does not fit the previously plotted KM estimates well!

# windows()
# tiff("Figure 13_10.tif", width = 4, height = 4, units = "in", res = 600, pointsize = 10, compression="lzw")

plot(xx, (1-dweibull(xx, theta1, omega)), type = "l", lwd = 2, col = "black",ylim=c(0,1),main="Weibull model estimate for Excellent (black) and Poor (grey) vitality classes",cex.main=0.5,xlab="Hours", ylab="Survival function")
lines(xx, (1-dweibull(xx, theta2, omega)), lwd = 2, col = "grey")

# dev.off()

### -------------------------------------------------------------------- ###


#!=========================================================================!#
#!====   Info Box 13.8 - Mixed distribution model - Worked example    =====!#
#!=========================================================================!#

# an example of model fit for a mixed distribution model (aka Cure-rate model) 
# using the optim function for a Weibul mixture model with vitality effects on theta and pi

## define the required variables
x=codKM$hours
y=codKM$censored
X=codKM$cond1 #dummy variable for vitality

## likelihood function used in the optimization
loglike=function(parms0, x, X, y){
  # Parse parameter vector:
  omega=exp(parms0[1])

  beta=parms0[2:5]
  
  linpred=beta[1] + beta[2]*X
  theta=exp(-linpred)
  
  # Linear predictor for mixture
  mu=beta[3] + beta[4]*X
  
  # Calculate survivor group probability from linear model:
  pi=1/(1+exp(-mu))
  
  # Calculate censored survival for each group:
  G_t=exp(-(theta*x)^omega)
  
  # Calculate density for each uncensored group:
  g=(omega * theta) * ((x*theta)^(omega-1)) * G_t
  
  # Calculate neg log-likelihood values:
  ll=-(((1-y)*(log(pi)+log(g))) + (y*log(1-pi + pi*G_t)))
  
  # Return sum of log-likelihood values:
  return(sum(ll))
}

## Define initial parameter vector:
parms0=c(log(1.0), 0, 0, 0,0) 

## Calculate test log-likelhood value for initial parameter vector to ensure the functions are defined
loglike(parms0, x, X, y)

## Perform minimization of neg log likelihood:
temp=optim(parms0, loglike, x = x, X = X, y = y,hessian=TRUE, control = list(fnscale = 1, trace = 3, maxit = 1500))
parms0=temp$par #parameters
parms=cbind(parms0, SE= sqrt(diag(solve(temp$hessian))) ) #parameters and standard errors

AIC=temp$value*2 + 2*nrow(parms)
AICc=AIC + (2*nrow(parms)*(nrow(parms)+1))/(length(y)-nrow(parms)-1)

### Plot results:

## 1.extract parameter estimates
xx=seq(0, 150, len = 150)
omega=exp(parms0[1])
X1=matrix(c(1,1,1,0),2,2,byrow=T)

linpred <- X1 %*% parms0[2:3]
theta=exp(-linpred)
mu <- X1 %*% parms0[4:5]
pi=1/(1+exp(-mu))

## calculate confidence intervals on pi
cov=solve(temp$hessian)
mu.se <- sqrt(diag(X1 %*% cov[4:5,4:5] %*% t(X1)))
pi.uci  <- exp(mu + 1.96 *mu.se) / (1 + exp(mu  + 1.96 *mu.se))
pi.lci  <- exp(mu - 1.96 *mu.se) / (1 + exp(mu  - 1.96 *mu.se))

## 2. plot the KM estimates to compare model fit
km_est <- survfit(Surv(hours, (1-censored)) ~ condition, data=codKM) #not that the function uses an indicator for survival rather than for censoring hence the use of (1-censored)

# windows()
# tiff("Figure 13_11.tif", width = 4, height = 4, units = "in", res = 600, pointsize = 10, compression="lzw")

plot(km_est, main="KM estimate (dash) and mixture model fit (solid) by vitality class",lty=2,cex=1,cex.main=0.5,
     xlab="Hours", ylab="Survival function", col= c("black","darkgrey"))
rm=matrix(nrow=2,ncol=2)
rm[1,1]=1
rm[1,2]=as.numeric(km_est$strata[1])
rm[2,1]=rm[1,2]+1
rm[2,2]=rm[1,2]+as.numeric(km_est$strata[2])

cond=1
km_est1=as.data.frame(cbind(time=km_est$time[rm[cond,1]:rm[cond,2]],upper=km_est$upper[rm[cond,1]:rm[cond,2]],lower=km_est$lower[rm[cond,1]:rm[cond,2]]))
with(km_est1,polygon(c(time, time[length(time):1]), c(lower, upper[length(time):1]), col =  rgb(0, 0, 0, alpha = 0.35), border = NA))

cond=2
km_est2=as.data.frame(cbind(time=km_est$time[rm[cond,1]:rm[cond,2]],upper=km_est$upper[rm[cond,1]:rm[cond,2]],lower=km_est$lower[rm[cond,1]:rm[cond,2]]))
with(km_est2,polygon(c(time, time[length(time):1]), c(lower, upper[length(time):1]), col =  rgb(0.5, 0.5, 0.5, alpha = 0.35), border = NA))


## 3. Plot curves for each injury code:

    # Weibull cumulative density function 
  dweibull=function(x, theta, omega){
    G=1 - exp(-(theta*x)^omega)
   return(G)
  }

points(xx, pi[1]*(1-dweibull(xx, theta[1], omega))+1-pi[1] , type = "l", lwd = 2, col = "black")
points(xx, pi[2]*(1-dweibull(xx, theta[2], omega))+1-pi[2] , type = "l", lwd = 2, col = "darkgrey")

# dev.off()


# tiff("cure rate 2.tif", width = 4.5, height = 4.5, units = "in", res = 600, pointsize = 10, compression="lzw")

plot(xx, pi[1]*(1-dweibull(xx, theta[1], omega))+1-pi[1] , type = "l", lwd = 2,ylim=c(0,1),xlab="Time",ylab='Survival probability')
points(xx, pi[2]*(1-dweibull(xx, theta[2], omega))+1-pi[2] , type = "l", lwd = 2, col = "red")

# dev.off()

### -------------------------------------------------------------------- ###


#!=========================================================================!#
#!==== Section 13.2.4.3	Generalized model for discard mortality data  =====!#
#!=========================================================================!#


# A generalized discard mortality model applied to data from Capizzano et al (2016; ICES J Mar Sci)
# Mixture-ditribution model with capture-handling mortality and natural mortality
# No effects of covariates


### Read in data
dat=read.csv("Capizzano et al 2016 cod.csv") #data from Capizzano et al (2016)

### Kaplan-Meier plot

# windows()
# tiff("Fig_all fish.tif", width = 4, height = 4, units = "in", res = 600, pointsize = 10, compression="lzw")

maxy=100
km_est <- survfit(Surv(event.time, (1-right.censor)) ~ 1, data=dat)
with(km_est,plot(time,surv,xlim=c(0,maxy),ylim=c(0.5,1),xaxs="i",yaxs="i", pch=1,lty=2,cex=1,xlab="Days", ylab="Survival probability", col= c("darkgrey")))
with(km_est,lines(time,surv,lty=1,col="darkgrey"))
with(km_est,polygon(c(time[time<=maxy], time[length(time[time<=maxy]):1]), c(lower[time<=maxy], upper[length(time[time<=maxy]):1]), col =  rgb(0, 0, 0, alpha = 0.1), border = NA))

#dev.off()

### subset the required data
x=dat$event.time
y=dat$right.censor
lcensor=dat$left.censor

### Define functions for model fitting

## Weibull density function
pweibull=function(x, theta, omega){
  G=exp(-(theta*x)^omega)
  g=(omega * theta) * ((x*theta)^(omega-1)) * G
  return(g)
}

## Weibull cumulative density function 
dweibull=function(x, theta, omega){
  G=1 - exp(-((theta*x)^omega))
  return(G)
}

## log-likelihood function for the generalized discard mortality model
loglike=function(parms0, x, y){
  # Parse parameter vector:
  theta=exp(-parms0[1])
  omega=exp(parms0[2])
  nu=parms0[3]
  mu=parms0[4]
  M=exp(parms0[5])
  
  # Calculate capture mort probability:
  tau=1/(1+exp(-nu))
  
  # Calculate survivor group probability from linear model:
  pi=1/(1+exp(-mu))
  
  # Calculate censored survival for each group:
  G1 = tau*(exp(-(M*x/365)-(theta*x)^omega))
  G2 = tau*exp(-M*x/365)
  
  # Calculate density for each uncensored group:
  g1 <- (((omega * theta) * ((x*theta)^(omega-1))) + (M/365)) * G1
  g2 <- (M/365) * G2
  
  # Calculate mixture density for uncensored group:
  h <- (pi * g1 + (1-pi) * g2)
  
  # Calculate mixture probability for right-censored group:
  H <- (pi * G1 + (1-pi) * G2)
  
  # Calculate mixture probability for left-censored group:
  H2 <- pi *(1-G1) + (1-pi) *(1-G2)
  
  # Calculate log-likelihood values:
  ll <- -(1-y)*(1-lcensor)*log(h) - y*log(H) - lcensor*log(H2)
  
  
  # Return sum of log-likelihood values:
  return(sum(ll))
}

### Model fitting routine

## Define initial parameter vector:
parms0=c(0.1, 0, 0, 0,log(0.2)) 

## Calculate test log-likelhood value for initial parameter vector:
loglike(parms0, x, y)

## Perform minimization of negative log-likelihood:
temp=optim(parms0, loglike, x = x, y = y,hessian=TRUE, control = list(fnscale = 1, trace = 3, maxit = 1500))
parms0=temp$par
parms=cbind(parms0, SE= sqrt(diag(solve(temp$hessian))) )
AIC=temp$value*2 + 2*nrow(parms)
(AICc=AIC + (2*nrow(parms)*(nrow(parms)+1))/(length(y)-nrow(parms)-1))

## estimate survival and 95%CI
cov=solve(temp$hessian)
library(MASS)
parmrand=mvrnorm(1000,mu=parms0,cov)

## Natural mortality
(M=exp(parms0[5]))
(M_CI=quantile(exp(parmrand[,5]),c(0.025,0.5,0.975)))

## capture and handling mortality (CH), 1-tau
(CHmort=1-(1/(1+exp(-(parms0[3])))))
(CI_CHmort=quantile(1-(1/(1+exp(-(parmrand[,3])))),c(0.025,0.5,0.975)))

## conditional post-release mortality (R), tau*pi
(Rmort=(1/(1+exp(-(parms0[3]))))*(1/(1+exp(-(parms0[4])))))
(CI_Rmort=quantile((1/(1+exp(-(parmrand[,3]))))*(1/(1+exp(-(parmrand[,4])))),c(0.025,0.5,0.975)))

## total CHR mortality
(total=1-1/(1+exp(-(parms0[3])))+(1/(1+exp(-(parms0[4]))))*(1/(1+exp(-(parms0[3])))))
(CI_total=quantile(1- (1/(1+exp(-(parmrand[,3])))) +(1/(1+exp(-(parmrand[,3]))))*(1/(1+exp(-(parmrand[,4])))),c(0.025,0.5,0.975)))

### Plot model fit against the Kaplan-Meier fit

## define plotting parameters
xmax=100
xx=seq(0, xmax, len = xmax+1)
theta=exp(-parms0[1])
omega=exp(parms0[2])
tau=1/(1+exp(-(parms0[3])))
pi=1/(1+exp(-(parms0[4])))

## Plot

# windows()
# tiff("capiz.tif", width = 4.5, height = 4.5, units = "in", res = 600, pointsize = 10, compression="lzw")

with(km_est,plot(time,surv,xlim=c(0,maxy),ylim=c(0.5,1),xaxs="i",yaxs="i", pch=1,lty=2,cex=1,xlab="Days", ylab="Survival probability", col= c("darkgrey")))
with(km_est,lines(time,surv,lty=1,col="darkgrey"))
with(km_est,polygon(c(time[time<=maxy], time[length(time[time<=maxy]):1]), c(lower[time<=maxy], upper[length(time[time<=maxy]):1]), col =  rgb(0, 0, 0, alpha = 0.1), border = NA))
points(xx, tau*(pi*(1-dweibull(xx, theta, omega))+1-pi)*(exp(-M*xx/365)) , type = "l", lwd = 2, col = "black")
legend(10,0.65,c("Kaplan-Meier","Model fit (max likelihood estimate)"),bt="n",cex=0.8,col=c("grey","black"),lty=c(1,1),lwd=3) 

# dev.off()

### -------------------------------------------------------------------- ###
### -------------------------------------------------------------------- ###
