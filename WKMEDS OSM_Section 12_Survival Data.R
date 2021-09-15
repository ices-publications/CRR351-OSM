#############################################################################
##### ICES Workshop on Methods for Estimating Discard Survival (WKMEDS) #####
#####  Online Supporting Materials for the Cooperative Research Report  #####
#####                                                                   #####
#####     Section 12 - Survival Data: format, structure & analysis      #####
#####                                                                   #####
#####                 Version 2.0      10th June 2021                   #####
#############################################################################

# Original code by M Breen & H Benoît 2017-02-10
# revised 2017-03-12, 2017-02-15, 2017-02-28, 2018-01-20, 2020-05-09, 2021-05-14
# Edited by M Breen 2021-06-10

# This code supports the examples and illustrations presented in section 12 of 
# the WKMEDS Cooperative Research Report (CRR).

### -------------------------------------------------------------------- ###


### Prep & Required Packages

rm(list=ls())

library("binom")
library("survival")



## Set Working Directory
# This code is written assumming the user has it and  the supporting data in a working directory.

setwd("E:/IMR/UM/ICES WKMEDS/CRR Report/WKMEDS CRR Final Edit/Online Support Materials/Chapter_12_Survival Data")

### -------------------------------------------------------------------- ###


#!=========================================================================!#
#!=======                   Cross Sectional Data                     ======!#
#!=======                       Section 12.2                         ======!#
#!=========================================================================!#

#!=========================================================================!#
#!======  12.2.1 The binomial distribution                           ======!#
#!=========================================================================!#


#===> Getting to Know the The Binomial Distribution <===#

library("binom")

?dbinom()  # for Binomial distribution help file

n <- 300
p <- 0.05

Value <- rep(0:n,1) # sets up a reference vector with integer values for X from 0 to n

## Generate & plot the density (cumulative probability) distribution
Prob <- dbinom(Value,n,p)  
plot(Dist~Value, type="l")

## Generate & plot the discrete probability distribution
Dist <- pbinom(Value, n, p) 
plot(Prob~Value, type="l")


### Binomial Probability Distribution Mean 

# Mean of a discrete probability function is µ = E[X] 
# i.e. the sum of the frequency (probability) times the value

Expected <- Value * Prob
Mean1 <- sum(Expected) 
Mean1

Mean2 <- n*p
Mean2

### Binomial Probability Distribution Variance
# Variance of a discrete probability function is E[(X-µ)^2] 
# i.e. the sum of the frequency (probability) times the squares of diff between the Value and mean

Var1 <- sum(Prob*(ifelse(Value<=n, Value - Mean1, 0))^2) 
Var1

# NB: the ifelse function is used to ensure only relevant values are used i.e. when value <= n, 
# otherwise Var will be artificially inflated because the vector structure will include values > n

Var2 <- n*p*(1-p)
Var2

SD <- sqrt(Var2)
SE <- SD/sqrt(n)



### Generating quantile values from binomial distribution

# These are synonimous with z values from Normal distribution & t values from Student's T distribution
# But they cannot be used to generate confidence intervals for the proportion statistic 
qbinom(0.975, n, p)
qbinom(0.025, n, p, lower.tail = TRUE)

### -------------------------------------------------------------------- ###


#!=========================================================================!#
#!=======    Plotting Binomial Discrete Probability Distribution     ======!#
#!=======                        Fig 12.02                           ======!#
#!=========================================================================!#

### Set up input values

n <- 50
p <- 0.2

Value <- rep(0:n,1) # sets up a reference vector with integer values for X from 0 to n

### Generate probabilty distribution
Prob <- dbinom(Value,n,p)  


### Bar Plot of Discrete Binomial Distribution with corresponding Normal Distribution

## Define "x" dimensions for Normal distribution
x_bar <- barplot(Prob)
df_x <- cbind(Value, x_bar)
df_x

## Define normal distribution based on binomial estimates for mean and sd
yv <- dnorm(Value, mean=n*p, sd=sqrt(n*p*(1-p)) ) 

## Set up Tiff file export
# to working directory
# tiff("Fig_12_02.tiff", units="cm", width=15, height=13, res=500)

# create positions for tick marks, one more than number of bars
at_tick <- seq_len(length(Value) + 1)

barplot(Prob, space = 0, axes = FALSE,
        ylim = c(0,max(Prob)+0.05)) 

# add y-axis
axis(side = 2, pos = -0.2)

# add x-axis with centered position, with labels, but without ticks.
axis(side = 1, at = seq_along(Value) - 0.5, pos = 0.01, tick = FALSE, labels = Value)

# Add axis titles
mtext("number of survivors", side = 1, line = 1.5)
mtext("probability", side = 2, line = 1.5)

## Smooth line 
# need to offset mean by 0.5 to fit
curve( (dnorm(x, mean=n*p+0.5, sd=sqrt(n*p*(1-p)) ) ),       
       col = "blue", lwd = 2, add=TRUE) 

##Close Tiff file export
# dev.off()

### -------------------------------------------------------------------- ###


#!=========================================================================!#
#!=======                        Fig 12.03                           ======!#
#!=======       Plotting the Relationship between se and p & n       ======!#
#!=========================================================================!#


## Set up baseline probability vectors

p_hat <- seq(0, 1, by = 0.001)
n_05 <- sqrt(p_hat*(1-p_hat)/5)
n_10 <- sqrt(p_hat*(1-p_hat)/10)
n_20 <- sqrt(p_hat*(1-p_hat)/20)
n_50 <- sqrt(p_hat*(1-p_hat)/50)


## Set up Tiff file export
# to working directory
# tiff("Fig_12_03.tiff", units="cm", width=15, height=13, res=500)

plot(NULL, xlim=c(0, 1), ylim=c(0, 0.25), frame.plot = FALSE, ylab="standard error", xlab=expression( hat(p) ) )
lines(p_hat, n_05, lty=1, lwd=1.5, col="blue")
lines(p_hat, n_10, lty=2, lwd=1.5, col="red")
lines(p_hat, n_20, lty=3, lwd=1.5, col="dark green")
lines(p_hat, n_50, lty=5, lwd=1.5, col="black")

### Add legend
## Legend text
Legend_txt <- c("n =  5", "n = 10", "n = 20", "n = 50")
Line_type <- c(1, 2, 3, 5)
Line_col <- c("blue", "red", "dark green", "black")
Line_width <- c(1.5, 1.5, 1.5, 1.5)
Point_type <- c(NA, NA, NA, NA)
legend("topright", Legend_txt, 
       cex = 0.9,
       pt.cex = 0,
       bty = "n", 
       col = Line_col, 
       lty = Line_type,
       lwd = Line_width,
       pch = Point_type)

##Close Tiff file export
# dev.off()

### -------------------------------------------------------------------- ###


#!=========================================================================!#
#!======  12.2.2 Confidence Intervals: a measure of uncertainty.     ======!#
#!=========================================================================!#

#!=========================================================================!#
#!======  Fig 12.05a - One Tailed Confidence Interval                ======!#
#!=========================================================================!#


### Plot Normal curve

## Set up input values
sigma <- 1
mu    <- 0

x  <- seq(-4, 4, by = 0.1)
y  <- ( 1/(sigma * sqrt(2*pi)) ) * ( exp(1)^( (-1 * ((x - mu)^2)) / (2*(sigma^2)) ) )

## Plot

## Set up Tiff file export
# to working directory
# tiff("Fig_12_05a.tiff", units="cm", width=15, height=13, res=500)

plot(x,y,type="l", lwd=2, col="blue", axes=F,  # turn axes off to align better - see below
     cex.lab = 1,
     ylim = c(0, 0.4),
     mgp = c(2,1,0),
     xlab = "Standard deviations from mean",
     ylab = "Probability Density") 
axis(1, pos=0, cex.axis=0.8)
axis(2, pos=-4.1, cex.axis=0.8)
segments(x0 = 0, y0 = 0, x1=0, y1=0.4, lty = 2, col= "grey")


### Plot one-sided 95% quantile

## Set up input values
Alpha <- 0.05
z_value_Alpha <- qnorm(1-Alpha)  # = 1.644854

lower.x <-  -4
upper.x <-  -z_value_Alpha

x=seq(lower.x, upper.x, by = 0.001) # Note need higher resolution than "by=0.1" to avoid slope in vertical
y  <- ( 1/(sigma * sqrt(2*pi)) ) * ( exp(1)^( (-1 * ((x - mu)^2)) / (2*(sigma^2)) ) )

## Plot
polygon(c(lower.x,x,upper.x), c(0,y,0), col="gray")
lines(x, y, col="blue", lwd=2)

### Annotate plot
text("a)", x=-5.3, y=0.45)   # now out of plotting range?!
text(expression( underline("One Sided 95% Quantile:") ), x=-2.3, y=0.36, cex=0.8)
text("-qnorm(1 - 0.05) = -1.64", x=-2.4, y=0.32, cex=0.8 )
text("95%", x=0, y=0.2, cex=1)
text("5%", x=-3, y=0.1, cex=1)

arrows(x0=-2, y0=0.03, x1=-3, y1=0.06, code=1, lwd=2, length=0.1, col="black", angle=20)
arrows(x0=-1.64, y0=0.1, x1=-2.2, y1=0.3, code=1, lwd=2, length=0.1, col="black", angle=20)

## Close Tiff file export
# dev.off()

### -------------------------------------------------------------------- ###


#!=========================================================================!#
#!======  Fig 12.05b - Two Tailed Confidence Interval                ======!#
#!=========================================================================!#


### Plot Normal curve

## Set up input values
sigma <- 1
mu    <- 0

x  <- seq(-4, 4, by = 0.1)
y  <- ( 1/(sigma * sqrt(2*pi)) ) * ( exp(1)^( (-1 * ((x - mu)^2)) / (2*(sigma^2)) ) )

## Plot

## Set up Tiff file export
# to working directory
# tiff("Fig_12_05b.tiff", units="cm", width=15, height=13, res=500)
plot(x,y,type="l", lwd=2, col="blue", axes=F,  # turn axes off to align better - see below
     cex.lab = 1,
     ylim = c(0, 0.4),
     mgp = c(2,1,0),
     xlab = "Standard deviations from mean",
     ylab = "Probability Density") 
axis(1, pos=0, cex.axis=0.8)
axis(2, pos=-4.1, cex.axis=0.8)
segments(x0 = 0, y0 = 0, x1=0, y1=0.4, lty = 2, col= "grey")


### Plot two-sided 95% quantile

## Set up input values
Alpha <- 0.05
z_value_Half.Alpha <- qnorm(1-Alpha/2)  # = 1.959964

lower.x1 <-  -4
upper.x1 <-  -z_value_Half.Alpha

x1=seq(lower.x1, upper.x1, by = 0.001) # Note need higher resolution than "by=0.1" to avoid slope in vertical
y1  <- ( 1/(sigma * sqrt(2*pi)) ) * ( exp(1)^( (-1 * ((x1 - mu)^2)) / (2*(sigma^2)) ) )

lower.x2 <-  z_value_Half.Alpha
upper.x2 <-  4

x2=seq(lower.x2, upper.x2, by = 0.1)
y2  <- ( 1/(sigma * sqrt(2*pi)) ) * ( exp(1)^( (-1 * ((x2 - mu)^2)) / (2*(sigma^2)) ) )

## Plot
polygon(c(lower.x1,x1,upper.x1), c(0,y1,0), col="gray")
lines(x1, y1, col="blue", lwd=2)

polygon(c(lower.x2,x2,upper.x2), c(0,y2,0), col="gray")
lines(x2, y2, col="blue", lwd=2)


### Annotate plot
text("b)", x=-5.3, y=0.45) # now out of plotting range?!
text("95%", x=0, y=0.2, cex=1)

text(expression( underline("Two Sided 95% Quantile:") ), x=-2.3, y=0.36, cex=0.8)
text("-qnorm(1 - 0.05/2) = -1.96", x=-2.4, y=0.32, cex=0.7 )
text("2.5%", x=-3, y=0.1, cex=1)

arrows(x0=-2.2, y0=0.02, x1=-3, y1=0.06, code=1, lwd=2, length=0.1, col="black", angle=20)
arrows(x0=-1.96, y0=0.06, x1=-2.2, y1=0.3, code=1, lwd=2, length=0.1, col="black", angle=20)

text(expression( underline("Two Sided 95% Quantile:") ), x=2.3, y=0.36, cex=0.8)
text("qnorm(1 - 0.05/2) = 1.96", x=2.4, y=0.32, cex=0.7 )
text("2.5%", x=3, y=0.1, cex=1)

arrows(x0=2.2, y0=0.02, x1=3, y1=0.06, code=1, lwd=2, length=0.1, col="black", angle=20)
arrows(x0=1.96, y0=0.06, x1=2.2, y1=0.3, code=1, lwd=2, length=0.1, col="black", angle=20)

## Close Tiff file export
# dev.off()

### -------------------------------------------------------------------- ###


#!=========================================================================!#
#!======  Table 12.1 - Different Binomal CIs                         ======!#
#!=========================================================================!#


binom.confint(x=70, n=100, 0.95, method = "all")



### -------------------------------------------------------------------- ###


#!=========================================================================!#
#!======  12.3.1 Survival/failure analysis: an introduction          ======!#
#!=========================================================================!#

#!=========================================================================!#
#!=======                        Fig 12.08                           ======!#
#!=======      Example of a Kaplan-Meier (KM) survival function      ======!#
#!=======          with a 95% confidence band (dashed lines)         ======!#
#!=========================================================================!#

# Figure 12.8.  Example of a Kaplan-Meier (KM) survival function, 
# with a 95% confidence band (dashed lines) using the data in tables 12.2 and 12.3.

### Load Data

Long_data_eg <- read.delim("E:/IMR/UM/ICES WKMEDS/CRR Report/WKMEDS CRR Final Edit/Online Support Materials/Chapter_12_Survival Data/Long_data_eg.txt")
View(Long_data_eg)
str(Long_data_eg)

attach(Long_data_eg)

### Generate KM Survival Function

library(survival)

## Surv() creates a survival object from the vectors "t" and "Censored"
# note - the surv() function uses an indicator for survival status, where dead = 1 and alive = 0.

my.surv <- Surv(t, Status) 

# alternatively, in this data, survival status can also be inferred from the censoring status,
# because all surviving individuals were right censored (i.e. Censored = 1).

# my.surv <- Surv(t, (1-Censored))  

# See section 13.2 where this coding format is used in worked examples 
# (i.e. Information Boxes 13.5, 13.6 and 13.7)


## survfit() fits a KM Survival Function to the data

my.fit <- survfit(my.surv ~ 1) 


### Plot the KM Survival Function with 95% confidence bands

## Set up Tiff file export
# to working directory
# tiff("Fig_12_08.tiff", units="cm", width=15, height=11, res=500)

plot(my.fit, main="KM Survival Function with 95% confidence band",
     xlab="time", ylab="survival function",
     conf.int = TRUE,
     bty = "l")

## Close Tiff file export
# dev.off()

## detach dataframe
detach(Long_data_eg)


### -------------------------------------------------------------------- ###


#!=========================================================================!#
#!======               12.4 Summarising survival data                ======!#
#!=========================================================================!#

#!=========================================================================!#
#!====== Fig 12_09 Importance of presenting effect size and CIs      ======!#
#!=========================================================================!#


## Define Wilson CIs

CI_A <- binom.wilson(x = 60, n = 70, 0.95)
CI_B <- binom.wilson(x = 49, n = 70, 0.95)
CI_C <- binom.wilson(x = 70, n = 100, 0.95)
CI_D <- binom.wilson(x = 975, n = 1500, 0.95)


## Combine into dataframe
str(CI_A)

DF_CIS <- rbind(CI_A, CI_B, CI_C, CI_D)
row.names(DF_CIS) <- c("A", "B","C", "D")
DF_CIS$Names <- c("A", "B","C", "D")
DF_CIS$line <- c(4, 3, 2, 1)

str(DF_CIS)

### Plot

## Plot without axes
plot(DF_CIS$mean, DF_CIS$line, axes=FALSE,
     ylab="",
     ylim= c(0.5, 4.5),
     xlab= "Survival",
     xlim= c(0,1),
     pch = 16,
     cex = 2,
     
     cex.lab=1.5,
     col.lab="grey")

### annotate plot

## with axes
axis(1, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), col.ticks="grey", col="black", col.axis="grey" )
axis(2, at=c(1, 2, 3, 4), labels=rev(DF_CIS$Names), pos = 0, las=2, tck=0, col="white" )
abline(v=0, col="black")

## CI limits
points(DF_CIS$lower , DF_CIS$line, pch="|")
points(DF_CIS$upper , DF_CIS$line, pch="|")
segments(x0= DF_CIS$lower, y0= DF_CIS$line, x1= DF_CIS$upper, y1= DF_CIS$line)

## Arbitrary Limit
abline(v=0.6, lty= 2, lwd=2, col="grey")


### Define Significance relative to nominal arbitary limit (H = 0.6)

Sig_A <- binom.test(60, 70, p=0.6, alternative = "greater")
Sig_B <- binom.test(49, 70, p=0.6, alternative = "greater")
Sig_C <- binom.test(70, 100, p=0.6, alternative = "greater")
Sig_D <- binom.test(975, 1500, p=0.6, alternative = "greater")




### -------------------------------------------------------------------- ###


#!=========================================================================!#
#!=======                        Fig 12.10                           ======!#
#!=======  Relationship between three different measures of effect:  ======!#
#!=======       risk difference, relative risk, and odds ratio       ======!#
#!=========================================================================!#


## Set up baseline probability vectors

p_1 <- seq(0.1, 1, by = 0.0001)
p_2 <- seq(0, 0.9, by = 0.0001)

RD <- p_1 - p_2
RR <- (1-p_2) / (1-p_1)
OR <- ( (1-p_2) / p_2 ) / ( (1-p_1) / p_1 )

## Set plotting limits <14
RR_lim <- ifelse(RR <= 14, RR, NA)
OR_lim <- ifelse(OR <= 14, OR, NA)


## Set up Tiff file export
# to working directory
# tiff("Fig_12_10.tiff", units="cm", width=15, height=13, res=500)

plot(NULL, xlim=c(0, 1), ylim=c(0, 14), frame.plot = FALSE, ylab="measure of effect", xlab=expression( hat(p)[1] ) )
lines(p_1, OR_lim, lty=1, lwd = 2, col="blue")
lines(p_1, RR_lim, lty=2, lwd = 2, col="red")
lines(p_1, RD, lty=5, lwd = 2, col="black")

### Add legend
## Legend text
Legend_txt <- c("Odds Ratio", "Relative Risk", "Risk Difference")
Line_type <- c(1, 2, 5)
Line_col <- c("blue", "red", "black")
Line_wid <- c(2, 2, 2)
Point_type <- c(NA, NA, NA)
legend("top", Legend_txt, 
       cex = 1.1,
       pt.cex = 0,
       bty = "n", 
       col = Line_col, 
       lty = Line_type,
       lwd = Line_wid,
       pch = Point_type)

##Close Tiff file export
# dev.off()


### -------------------------------------------------------------------- ###
### -------------------------------------------------------------------- ###

