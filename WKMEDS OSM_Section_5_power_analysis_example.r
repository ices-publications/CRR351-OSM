#############################################################################
##### ICES Workshop on Methods for Estimating Discard Survival (WKMEDS) #####
#####  Online Supporting Materials for the Cooperative Research Report  #####
#####                                                                   #####
#####        Section 5.3 - Experimental Design: Power Analysis          #####
#####                                                                   #####
#####                 Version 2.0      14th May 2021                    #####
#############################################################################

# Contributors: Chun Chen, Bob van Marlen and Mike Breen

# This code supports the power analysis examples presented in section 5.3 of 
# the WKMEDS Cooperative Research Report (CRR).


#!=========================================================================!#
#!===============     Power Analysis: Analytical Approach     =============!#
#!===============             Section 5.3.1                   =============!#
#!=========================================================================!#

# Code by Mike Breen (Institute of Marine Research (IMR), Norway).


### Required packages

library(binom)

### Example - estimating minumim sample size

binom.power(p.alt = 0.85, n = 50, p = 0.75, 
            alpha = 0.05,                           # set sig.level to 5%
            phi = 1,                                # no over-dispersion
            alternative = "greater",                
            method = "asymp")

### -------------------------------------------------------------------- ###

### Code to generate Figure 5.2

## Input parameters

  alpha <- 0.05                                     # set sig.level to 5%
  phi <- 1                                          # no over-dispersion
  p <- 0.75
  p.alt <- 0.85

  n = 1:200                                         # Range of sample sizes
  
#Set-up receiving vector
power.out <- rep("NA", length(n))

## Loop to generate statistical power estimates from a range of sample sizes
for (i in 1:length(n)) {
    power <- binom.power(p.alt = p.alt, n = i, p = p, alpha = alpha, phi = phi, alternative = "greater", method = "asymp")
    power.out[i] <- as.numeric(power[1])
    }

## Identify minimum required sample size
power.out <- as.numeric(power.out)
power.plot <- data.frame(n, power.out)
power.OK <- (power.plot[power.plot$power.out >= 0.8 , ])
sample.min <- as.numeric(min(power.OK$n))

## Plotting vectors for minimum sample size
x <- c(-5, sample.min, sample.min)
y <- c(0.8, 0.8, 0)

## Plot Statistical Power vs Sample Size
plot(n, power.out,
     type="l", lty=1, lwd=2, col="black",
     xlab="Sample size (n)", mgp=c(2.25,1,0), cex.axis=1, cex.lab=1.2,
     ylab=expression(paste("Statistical Power (1 - ", beta, ")", sep="")) )
lines(x,y, type="l", lty=2, lwd=1, col="black")


setwd("E:/IMR/UM/ICES WKMEDS/CRR Report/WKMEDS CRR Final Edit/Online Support Materials/Chapter_05_Experimental Design")

tiff(filename="Fig_5.2.tiff", width=800, height=500)

plot(n, power.out,
     type="l", lty=1, lwd=2, col="black",
     xlab="Sample size (n)", mgp=c(2.25,1,0), cex.axis=1.5, cex.lab=2,
     ylab=expression(paste("Statistical Power (1 - ", beta, ")", sep="")) )
lines(x,y, type="l", lty=2, lwd=1, col="black")

dev.off()


### -------------------------------------------------------------------- ###


#!=========================================================================!#
#!================= Monte Carlo Simulation based Power Analysis ===========!#
#!=================             Section 5.3.2                   ===========!#
#!=========================================================================!#

# Original R-code by Chun Chen of WMR, IJmuiden, the Netherlands.
# Text translated into English by Bob van Marlen (WMR). 
# Edited by Mike Breen (IMR, Norway).


### Parameters:
# b         : The true (assumed) effect for treatment as compared to control
# n         : Number of fishes per survival tank (survival tank = experimental unit)
# J         : Number of replicates per treatment cell
# b_0       : Interpcet in the GLMM model
# sigma.a   : standard deviation of the random intercept in the GLMM model
# POWER     : statistical power of the study, e.g. 0.8
# SIG_LEVEL : statistical significance level in the test, e.g. 0.05  
# Nsim      : Number of simulations

### -------------------------------------------------------------------- ###


# Given the number of fishes (n) in each of the experiment setting for the storage factor(cage and tank)
# and the number of replicates per experiment setting, the power of the factor origin is calculated. 
# With a desired power (e.g. 0.8) one can also determine the
# number of replicates needed by plotting power vs. number of replicates.

# The values used are based on survival experiments in 2013, see:
  # van Marlen, B., Nijman, R., Groeneveld, K., Goudswaard, K., Uhlmann, S.S., 2013.
  # Praktijk Netwerk Discards Zuid - Verkennend onderzoek aan visies omtrent
  # discardvermindering en overleving van ondermaatse vis. IMARES Report C126/13. p. 53.  
# Also see Report Uhlmann et al., xx]

# The experiment consists of taking live fish from differing spots of the processing
# line, e.g from the deckbin and from the belt, and storing fish in two different ways:
# in underwater cages and in tanks on-board.

# As we have 1 factor (origin) with four levels.
# levels are: hopper, begindiscard, enddiscard, and control
# J is the number of replicates for each origin setting, and we wish to determine the
# size of J by this power analysis.

#!===========================================================================!#

# install.packages("lme4") # Not needed when already done
library(lme4)

## Simulated Data for power analysis
simulate_data <- function(J, n, hp) {
  
  ## J = Number of replicates in each experiment setting (hopper, begindiscard, enddiscard, and control
  ## n = Number of fishes in a experiment
  
  ##==experiment setting: hopper ======= reference = control
  
  y1 <- NA
  y2 <- NA
  y3 <- NA
  y4 <- NA
  y1 <- hp["mu.a"]+hp["b.hopper"]                                             ##---fixed effect from hopper---
  y2 <- rnorm(J, 0, hp['sigma.a'])                                            ##---random effect among the experiments---
  y3 <- y1+y2                                                                 ##---establish the linear function
  P  <- 1/(1+(exp(-y3)))                                                      ##---apply the logistic transform---
  y4 <- rbinom(J, n, P)                                                      ##---for each replicate, the outcome (e.g. number of survial) follows a binomial distribution binomial(n, P)---
  temp1  <- data.frame(n.alive=y4, n.dead=n-y4, origin=rep("hopper",J), expID=seq(1,J), Prob=P, fix_term=rep(y1,J), ran_term=y2)
  ##---expID gives the index of the replicate
  
  ##==experiment setting: begindiscard======= reference= control
  
  y1 <- NA
  y2 <- NA
  y3 <- NA
  y4 <- NA
  y1 <- hp["mu.a"]+hp["b.begindiscard"]   
  y2 <- rnorm(J, 0, hp['sigma.a'])                 
  y3 <- y1+y2
  P  <- 1/(1+(exp(-y3)))
  y4 <- rbinom(J, n, P)    
  temp2  <- data.frame(n.alive=y4, n.dead=n-y4, origin=rep("begindiscard",J), expID=seq((J+1),2*J), Prob=P, fix_term=rep(y1,J), ran_term=y2)
  
  ##==experiment setting: enddiscard======= reference= control
  
  y1 <- NA
  y2 <- NA
  y3 <- NA
  y4 <- NA
  y1 <- hp["mu.a"]+hp["b.enddiscard"]   
  y2 <- rnorm(J, 0, hp['sigma.a'])                 
  y3 <- y1+y2
  P  <- 1/(1+(exp(-y3)))
  y4 <- rbinom(J, n, P)    
  temp3  <- data.frame(n.alive=y4, n.dead=n-y4, origin=rep("enddiscard",J), expID=seq(2*J+1,3*J), Prob=P, fix_term=rep(y1,J), ran_term=y2)
  
  ##==experiment setting: control (reference) ======= 
  
  y1 <- NA
  y2 <- NA
  y3 <- NA
  y4 <- NA
  y1 <- hp["mu.a"]                                             
  y2 <- rnorm(J, 0, hp['sigma.a'])                                            
  y3 <- y1+y2                                                                
  P  <- 1/(1+(exp(-y3)))                                                      
  y4 <- rbinom(J, n, P)                                                        
  temp4  <- data.frame(n.alive=y4, n.dead=n-y4, origin=rep("control",J), expID=seq(3*J+1,4*J), Prob=P, fix_term=rep(y1,J), ran_term=y2)
  
  ydata <- rbind(temp1, temp2, temp3, temp4)
  
}   # End of function simulate_data



#!===========================================================================!#

##----define the expected effect to be detected:                                                                                                                                               

hp        <- c(log(0.8), log(0.5), log(0.2), 1, 0.1)
names(hp) <- c("b.hopper", "b.begindiscard", "b.enddiscard", "mu.a", "sigma.a") 

# This is what one wants to achieve.
# 
# hp[1] ("b.hopper") is the coefficient for the relative effect of hopper compared to control (control is the reference group in the model).
# hp[2] ("b.begindiscard") is the coefficent for the relative effect of begindiscard vs. control.
# hp[3] ("b.enddiscard") is the coefficent for the relative effect of enddiscard vs. control.

# A negative coefficient means a negative contribution to survival (if the outcome is survival) as compared to the control (reference) group, 
#     e.g. fishes coming from the hopper have a lower survival than from the control group. 

# hp[4] ("mu.a") is the intercept in the model, but this is not very relevant to the power of the origin factor. One can use other values instead. 
#     In reality, it is associated with the survival rate of the control group. mu.a=0 indicates no death from the control group, but in reality it is usually not true.
# hp[5] ("sigma.a") is the standard deviation of the random effect from fishes among the different experiments.

# So the message here is that one needs to have some a-priori knowledge of effects, for instance conduct a pilot study.
# Alternatively one can try different values of the coefficients and find the corresponding coefficient/number of replicates that is both biologically sound and logistically afforable.       

#!===========================================================================!#

## Define range of replicates per treatment
# J <- c(2, 5, 10, 20, 30, 40, 50, 60, 70, 80)

# One can choose other values for J here, in this case we take 2, 5, 10.
J <- c(2, 5)

power_hopper  <- rep(NA, length(J))
power_begindiscard  <- rep(NA, length(J))
power_enddiscard  <- rep(NA, length(J))
power_origin  <- rep(NA, length(J))
# power_hopper refers to the power of detecting a non-equal survival between hopper and control (reference group). Thus, given the coefficient hp[1] and hp[5] 
#     parameters in th data dsitribution, what's the probability of detecting a non-equal survival between hopper and control.
# Similarly, power_begindiscard or power_enddiscard refer to the power of detecting a non-equal survival between begindiscard or enddiscard and control, respectively.
# power_origin refers to the power of detecting a non-equal survival among all levels of origin. 

for (i in 1:length(J)) {
  Nsim <- 1000
  P_hopper <- rep(NA, Nsim)
  P_begindiscard <- rep(NA, Nsim)
  P_enddiscard <- rep(NA, Nsim)
  P_origin <- rep(NA, Nsim)
  
  for (k in 1:Nsim) {
    
    #-----simulate data----
    dat <- simulate_data(J[i], n=10, hp)                                          ##---n: number of fishes per experiment
    dat$origin <- factor(dat$origin, levels = c("control", "hopper", "begindiscard", "enddiscard"))  ##----change the order of the levels in factor origin, so that control is the reference level
    
    #-----fit the model----
    fit <- glmer(cbind(n.alive, n.dead) ~ origin + (1|expID), family=binomial, data=dat,
                 control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)) )
    
    #-----extract p-values of comparing each level of origin with the reference level control----
    vc <- vcov(fit, useScale = FALSE)
    b  <- fixef(fit)
    se <- sqrt(diag(vc))
    z  <- b / sqrt(diag(vc))
    
    P_hopper[k]  <- 2 * (1 - pnorm(abs(z["originhopper"])))
    P_begindiscard[k] <- 2 * (1 - pnorm(abs(z["originbegindiscard"])))
    P_enddiscard[k] <- 2 * (1 - pnorm(abs(z["originenddiscard"])))
    
    #-----fit the model without the origin factor----
    fit0 <- glmer(cbind(n.alive, n.dead) ~ (1|expID), family=binomial, data=dat,
                  control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)) )
    #-----extract p-values of comparing all levels of origin as a whole----
    P_origin[k] <- anova(fit, fit0)[,"Pr(>Chisq)"][2]
    
  } # end of k-loop
  
  ##---compute the power for hopper vs. control when the number of replicates is J[i], thus the probability of detecing a non-equal survival between hopper and control---
  ##---the significance level is 0.05---
  power_hopper[i]  <- sum(P_hopper<0.05, na.rm = TRUE)/length(P_hopper[!is.na(P_hopper)])
  
  ##---compute power for begindiscard vs. control when the number of replicates is J[i]
  power_begindiscard[i] <- sum(P_begindiscard<0.05, na.rm = TRUE)/length(P_begindiscard[!is.na(P_begindiscard)])
  
  ##---compute power for enddiscard vs. control when the number of replicates is J[i]
  power_enddiscard[i] <- sum(P_enddiscard<0.05, na.rm = TRUE)/length(P_enddiscard[!is.na(P_enddiscard)])
  
  ##---compute power for hopper vs. begindiscard vs. enddiscard vs. control when the number of replicates is J[i]
  power_origin[i] <- sum(P_origin<0.05, na.rm = TRUE)/length(P_origin[!is.na(P_origin)])
} ##---end of loop for i


### -------------------------------------------------------------------- ###


############################
### Save Output as table ###
############################

## Sigma 0.1

N_replicates <- N

Data_save <- data.frame(N_replicates, power_begindiscard, power_enddiscard, power_hopper, power_origin)

print(Data_save)

write.csv(Data_save, file = "Data_save_sigma_0_1.csv")



## Sigma 0.5

N_replicates <- N
power_begindiscard1 <- power_begindiscard
power_enddiscard1 <- power_enddiscard
power_hopper1 <- power_hopper
power_origin1 <- power_origin

Data_save <- data.frame(N_replicates, power_begindiscard1, power_enddiscard1, power_hopper1, power_origin1)

print(Data_save)

write.csv(Data_save, file = "Data_save_sigma_0_5.csv")



Data_save_sigma_0_5 <- read.csv("E:/IMR/UM/ICES WKMEDS/CRR Report/WKMEDS CRR Final Edit/Online Support Materials/Chapter_05_Experimental Design/Data_save_sigma_0_5.csv")
View(Data_save_sigma_0_5)


power_begindiscard1 <- Data_save_sigma_0_5$power_begindiscard1
power_enddiscard1 <- Data_save_sigma_0_5$power_enddiscard1
power_hopper1 <- Data_save_sigma_0_5$power_hopper1
power_origin1 <- Data_save_sigma_0_5$power_origin1


### -------------------------------------------------------------------- ###


########################
### plot the results ###
########################


setwd("E:/IMR/UM/ICES WKMEDS/CRR Report/WKMEDS CRR Final Edit/Online Support Materials/Chapter_05_Experimental Design")

tiff(filename="myplot01.tiff", width=800, height=800)

plot(J, power_hopper, type="o", lwd=2, ylim=c(0,1), col="black", pch=15, cex=1.5,
      main=expression(paste(sigma[a]," = 0.1, n = 10", sep="")), cex.main=2,
      xlab="Number of replicates (J)", mgp=c(2.5,1,0), cex.axis=1.5, cex.lab=2,
      ylab=expression(paste("Statistical Power (1 - ", beta, ")", sep="")) )
lines(J, power_begindiscard, type="o", lwd=2, col="black", pch=16, cex=1.5)
lines(J, power_enddiscard, type="o", lwd=2, col="black", pch=17, cex=1.5, lty=2)
lines(J, power_origin, type="o", lwd=2, col="black", pch=1, cex=1.5)
abline(h=0.8, lty=3)
legend(x=25, y=0.22, bty="n",pch=c(15, 16, 17, 1), cex=2,legend=c(expression(paste("hopper vs. control: ", exp(b), "=0.8", sep="")), 
                                                                  expression(paste("begindiscard vs. control: ", exp(b), "=0.5", sep="")), 
                                                                  expression(paste("enddiscard vs. control: ", exp(b), "=0.2", sep="")), 
                                                                  expression(paste("all: all ", b, sep=" "))), lwd=2, 
                                                                  lty=c(1,1,2,1)
)  

dev.off()


tiff(filename="myplot02.tiff", width=800, height=800)

plot(J, power_begindiscard, type="o", lwd=2, ylim=c(0,1), col="black", pch=15, cex=1.5,
     xlab="Number of replicates (j)",mgp=c(2.5,1,0), cex.axis=1.5, cex.lab=2,
     ylab=expression(paste("Statistical Power (1 - ", beta, ")", sep="")), 
     main=expression(paste("begindiscard vs. control: ", exp(b), " = 0.5,", "n = 10", sep="")), cex.main=2)
lines(J, power_begindiscard1, type="o", lty=2, lwd=2, col="black", pch=16, cex=1.5)
abline(h=0.8, lty=3)
legend(x=60, y=0.22, bty="n",
       pch=c(15, 16, 17, 1), 
       cex=2,lwd=2,lty=c(1,1,2,1),
       legend=c(expression(paste(sigma[a], " = 0.1", sep="")), 
                expression(paste(sigma[a], " = 0.5", sep="")))
)  

dev.off()



### -------------------------------------------------------------------- ###








