###########################################################################################################
###########################################################################################################
### Data supporting the short-term health effects of temperature and air pollution in Valencia, Spain.  ###
### Data in Brief, under review (19/04/2022)                                                            ###
###                                                                                                     ###
### Carmen Iniguez (carmen.iniguez@uv.es)                                                               ###
### Ferran Ballester (ballester_fer@gva.es)                                                             ###
### Aurelio Tobias(aurelio.tobias@idaea.csic.es)                                                        ###
###                                                                                                     ###
### valencia0107.R                                                                                      ###
### Last update: 27/04/2022                                                                             ###
###########################################################################################################
###########################################################################################################

# Remove objects.
rm(list = ls())

# Install packages.
install.packages(c("splines","tsModel","dlnm"))

# Load  packages.
library(splines)
library(tsModel)
library(dlnm)

###########################################################################################################
### Load Valencia dataset.                                                                              ###
###########################################################################################################

# Load data.
data <- read.csv('valencia0107.csv', sep=",")

# List first rows of time-series data.
head(data)

# Formatting date.
data$date <- as.Date(data$date)

# Generating time variables.
data$time <- seq(nrow(data))
data$dow  <- weekdays(data$date)

###########################################################################################################
### Figure 1.                                                                                           ###
### Daily counts of all-cause mortality and mean temperature (?C) in Valencia, Spain, 2001-2007.        ###
###########################################################################################################

# Plot time-series.
pdf("Figure 1.pdf", width=15, height=10)
par(mex=0.8,mfrow=c(2,1))

  # Deaths. 
  plot(data$date, data$all, 
      type="l", col= "blue",
      ylim=c(0,50),
      ylab="Num. of deaths", xlab="Date")

  # Temperature.
  plot(data$date, data$tmax, 
       type="l", col= "red",
       ylim=c(5,42),
       ylab="Temperature (ºC)", xlab="Date")

dev.off()
layout(1)

###########################################################################################################
### Figure 2.                                                                                           ###
### Relative risk (RR) of mortality along daily temperature and lag dimension with reference at 21ºC;   ###
### plot of RR for temperature of 28ºC at specific lags (bottom-left panel);                            ###
### and plot of RR at lag 0 of temperature (bottom-right panel).                                        ###
###########################################################################################################

# Natural cubic splines for time-trend and seasonality with 10 df/year.
numyears <- length(unique(data$year))
df <- 10
spl <- ns(data$time, df*numyears)

# Crossbasis for temperature.
pct <- quantile(data$tmean, prob=c(.10,.75,.90),na.rm=T)
varknot <- pct[c(1,2,3)] # knots at 10th, 75th, and 90th
klag <- logknots(21, 3)
cb.temp <- crossbasis(data$tmean, lag=21, 
                      argvar=list(fun="ns", knots=varknot),
                      arglag=list(fun="ns", knots=klag) )
summary(cb.temp)

# Fit quasi-Poisson regression model.
model <- glm(all ~ spl + cb.temp + factor(dow) + factor(hol),  data=data, family=quasipoisson)
pred <- crosspred(cb.temp, model, by=1) 

# Get prediction centered at the MMT.
mmt <- pred$predvar[which.min(pred$allRRfit)] 
mmt
predcen <- crosspred(cb.temp, model, cen=mmt, by=1) 

# Plot DLNM for temperature. 
pdf("Figure 2.pdf", width=15, height=20)
l <- layout(matrix(c(1, 1, 2, 3, 4, 5), 
                   nrow = 3,
                   ncol = 2,
                   byrow = TRUE), heights = c(10,5,5))

  # 3D surface plot.
  temp3d <- plot(predcen, shade=0.15, col=grey(1),
                 ylim=c(0,21), xlim=c(2,32), zlim=c(0.95,1.15),
                 ylab="Lag", xlab="Temperature (ºC)", zlab="RR")
  lines(trans3d(x= 7, y=0:21, z=predcen$matRRfit[as.character( 7),], pmat=temp3d), lwd=5, col="blue")
  lines(trans3d(x=28, y=0:21, z=predcen$matRRfit[as.character(28),], pmat=temp3d), lwd=5, col="blue")
  lines(trans3d(x=predcen$predvar, y=0, z=predcen$matRRfit[,"lag0"], pmat=temp3d), lwd=5, col="red")
  lines(trans3d(x=predcen$predvar, y=14, z=predcen$matRRfit[,"lag14"], pmat=temp3d), lwd=5, col="red")
  
  # Plot for cold at specific lags.
  plot(predcen, var=7, ylim=c(0.9,1.2), lwd=4,
       main="Temperature at 7 ºC", ylab="RR", col="blue")
  box(lty = 1)
  
  # Plot for lag 0 at each temperature.
  plot(predcen, lag=0, ylim=c(0.9,1.2), lwd=4,
       main="Lag 0", ylab="RR", xlab="Temperature (ºC)", col="red")
  box(lty = 1)

  # Plot for heat at specific lags.
  plot(predcen, var=28, ylim=c(0.9,1.2), lwd=4,
       main="Temperature at 28 ºC", ylab="RR", col="blue")
  box(lty = 1)
  
  # Plot for lag 12 at each temperature.
  plot(predcen, lag=14, ylim=c(0.9,1.2), lwd=4,
       main="Lag 14", ylab="RR", xlab="Temperature (ºC)", col="red")
  box(lty = 1)
  
dev.off()
layout(1)

###########################################################################################################
### Figure 3.                                                                                           ###
### Overall cumulative exposure-response association between daily temperature and mortality across     ###
### all lags, with related temperature distribution.                                                    ###
###########################################################################################################

# Get relative risk (RR) for cold at 2.5th percentile of temperature vs.MMT.
cold <- round(quantile(data$tmean, prob=c(.025), na.rm=T))
cold 
target <-as.character(cold)
RR     <- predcen$allRRfit[target] 
RRlow  <- predcen$allRRlow[target] 
RRhigh <- predcen$allRRhigh[target] 
cbind(RR, RRlow, RRhigh)

# Get relative risk (RR) for heat at 97.5th percentile of temperature vs MMT.
heat <- round(quantile(data$tmean, prob=c(.975), na.rm=T))
heat
target <-as.character(heat)
RR     <- predcen$allRRfit[target] 
RRlow  <- predcen$allRRlow[target] 
RRhigh <- predcen$allRRhigh[target] 
cbind(RR, RRlow, RRhigh)

# Plot overall cumulative exposure-response.
pdf("Figure 3.pdf", width=10, height=8)

  plot(predcen, "overall",
       ylim=c(0.8,3), xlim=c(2,32), 
       lwd=2,
       ylab="RR", xlab="Temperature (ºC)")
  abline(v=c(mmt, cold, heat), lty=(c(1,2,2)))
  box(lty = 1)

dev.off()

###########################################################################################################
###########################################################################################################
###                                                                                                     ###
###                                     End of script file                                              ###
###                                                                                                     ###
###########################################################################################################
###########################################################################################################
