###############################################################################
# This code, writtin by Benjamin Golas, will take output from Bayesian 
# parameterization and fit those
# models back into the Hayman/Haase models for comparison of survival
# estimates
###############################################################################

rm(list=ls())

library(tidyverse)
library(reshape2)
library(ghibli)
library(cowplot)


setwd("C:/Users/bengo/Documents/WNS/BATPACK_ANALYSIS/")
df <- as.data.frame(read.csv2("EnsembleGuidedUmbrellaIndividual.csv"))[,-1]

T.a.values <- rep(seq(-2,20,0.25), length(seq(0.4,1,0.025)))
RH.values <- rep(seq(0.4,1,0.025), each=length(seq(-2,20,0.25)))

#testing purposes
# T.a <- -2
# RH <- 1




T.HIB <- function(T.a,RH,df){

  
  
df <- as.data.frame(df) 
mu.T.a <- 6.005687 #scale from data
sd.T.a <- 4.135221
T.a.scale <- (T.a-mu.T.a)/sd.T.a
  
#Bat parameters
C.t <- df$mu.C.t
T.tor.min <- df$mu.T.tor.min
TMR.min <- df$mu.TMR.min
rEWL.body <- df$mu.rEWL.body
t.tor.max <- df$mu.t.tor.max

#hierarchical parameters
b1 <- df$beta.1.
b2 <- df$beta.2.
b3 <- df$beta.3.
b4 <- df$beta.4.
g1 <- df$g.1.
g2 <- df$g.1.
g3 <- df$g.1.
g4 <- df$g.1.

C <- 0.2
percent.fat <- 0.3
S <- 0.131
T.eu <- 35
WR <- 90
M.body <- 20
C.eu <- 0.20
RMR <- 1.13
T.lc <- 31.35
t.eu <- 1 #median from all of our bats




rEWL.wing <- rEWL.body*0.33/0.1

M.fat <- M.body*percent.fat
M.lean <- M.body - M.fat #total body mass

SA.body <- 10*M.body^(2/3) #body surface area
SA.wing <- SA.body*19.68/39.36 #wing surface area, per Haase2019

T.tor <- pmax(T.a,T.tor.min)
Q10 <- 1.6+0.26*T.a-0.006*T.a^2

t.cool <- (log(T.eu-T.tor)/log(10))/(C*M.body^0.67*(log(T.eu-T.a)/log(10))/(S*M.body))
Q10.cool <- 3.82-0.507*log10(M.body)
E.cool <- t.cool*(TMR.min +(RMR*Q10.cool^((T.a-T.eu)/10)))
t.warm <- (T.eu-T.a)/WR
E.warm <- S*(T.eu-T.a)+t.warm*(C.eu*(T.eu-T.a))

E.eu <- t.eu*(RMR+C.eu*(T.lc-T.a))

#Calculate d.WVP
WVP.bat <- .611*exp((17.503*T.tor)/(T.tor+240.97))
WVP.a <- RH*(.611*exp((17.503*T.a)/(T.a+240.97)))
d.WVP <- max(0.0001, WVP.bat-WVP.a)

#Hayman model
t.tor.mark <- ifelse(rep(T.a,length(T.tor.min))<=T.tor.min, 1,0)
rho <- exp(g1+g2*T.a.scale+g3*d.WVP+g4*T.a.scale*d.WVP)/(1+exp(g1+g2*T.a.scale+g3*d.WVP+g4*T.a.scale*d.WVP))
t.tor <- rho*(t.tor.mark*(t.tor.max/(1+(T.tor.min-T.a)*(C.t/TMR.min))) + (1-t.tor.mark)*(t.tor.max/Q10^((T.a-T.tor.min)/10)))
E.tor <- t.tor.mark*(TMR.min +(T.tor.min-T.a)*C.t) + (1-t.tor.mark)*(TMR.min*Q10^((T.a-T.tor.min)/10))

E.bout <- E.cool + E.tor*t.tor + E.warm + E.eu
E.hib <- M.fat*39.3*100/20.1

t.hib <- (E.hib/E.bout)*(t.tor+t.cool+t.warm+t.eu)/24/30

#Now with disease
t.winter <- 5000
T.Pd.min <- 0
T.Pd.max <- 19.7
B.Pd.1 <- 1.15*10^(-3)
B.Pd.2 <- 0.27
mu.Pd.1 <- 1.51*10^(-4)
mu.Pd.2 <- -9.92*10^(-3)

Pd <- max((B.Pd.1*(T.tor-T.Pd.min)*(1-exp(B.Pd.2*(T.tor-T.Pd.max)))*((mu.Pd.1*RH*100)/(1+mu.Pd.2*RH*100))*t.winter),1)
t.tor.Pd <- rho*(t.tor.mark*(t.tor.max/(1+(T.tor.min-T.a)*(C.t/TMR.min))) + (1-t.tor.mark)*(t.tor.max/Q10^((T.a-T.tor.min)/10)))/Pd

E.bout.Pd <- E.cool + E.tor*t.tor.Pd + E.warm + E.eu

t.hib.Pd <- (E.hib/E.bout.Pd)*(t.tor.Pd+t.cool+t.warm+t.eu)/24/30

output <- c(T.a,RH,mean(t.hib),quantile(t.hib,probs=c(0.025,0.5,0.975)),mean(t.hib.Pd),quantile(t.hib.Pd,probs=c(0.025,0.5,0.975)))
# output <- c(T.a,RH,mean(t.hib),quantile(t.hib,probs=c(0.025,0.5,0.975)))
return(output)

}#end function



library(doParallel)
detectCores()
c <- makeCluster(11)
registerDoParallel(c)
getDoParWorkers()

results <- foreach(i = 1:length(T.a.values), .combine= cbind) %dopar% {
  T.a <- T.a.values[i]
  RH <- RH.values[i]
  ind.result <- T.HIB(T.a=T.a,RH,df=df)
  return(ind.result)
}

df.bat1 <- as.data.frame(t(results))
colnames(df.bat1) <- c("Temp","RH","Mean","CILow","Median","CIHigh","Mean.Pd","CILow.Pd","Median.Pd","CIHigh.Pd")
# colnames(df.bat1) <- c("Temp","RH","Mean","CILow","Median","CIHigh")


df.data <- as.data.frame(read.csv("DataAverages.csv")[,-1])
df.data$bout.length <- df.data$bout.length/24
df.data.RH <- df.data$WVP.a/(.611*exp((17.503*df.data$T.a)/(df.data$T.a+240.97)))
df.data <- cbind(df.data,round(df.data$T.a*4,)/4, round(df.data.RH*40,)/40)
colnames(df.data)[c(5,6)] <- c("T.a.round","RH.round")
df.bat1$RH <- round(df.bat1$RH,3)
df.data$RH.round <- round(df.data$RH.round,3)
# df.bat1$Temp <- round(df.bat1$Temp,2)
survive.low <- rep(NA, length(df.data[,1]))
survive.mean <- rep(NA, length(df.data[,1]))
survive.high <- rep(NA, length(df.data[,1]))
survive.low.Pd <- rep(NA, length(df.data[,1]))
survive.mean.Pd <- rep(NA, length(df.data[,1]))
survive.high.Pd <- rep(NA, length(df.data[,1]))
for(i in 1:length(df.data[,1])){
  survive.low[i] <- ifelse(df.bat1$CILow[which(df.bat1$Temp==df.data$T.a.round[i] & df.bat1$RH==df.data$RH.round[i])]>=4780/24/30, "survive","die")
  survive.mean[i] <- ifelse(df.bat1$Mean[which(df.bat1$Temp==df.data$T.a.round[i] & df.bat1$RH==df.data$RH.round[i])]>=4780/24/30, "survive","die")
  survive.high[i] <- ifelse(df.bat1$CIHigh[which(df.bat1$Temp==df.data$T.a.round[i] & df.bat1$RH==df.data$RH.round[i])]>=4780/24/30, "survive","die")
  survive.low.Pd[i] <- ifelse(df.bat1$CILow.Pd[which(df.bat1$Temp==df.data$T.a.round[i] & df.bat1$RH==df.data$RH.round[i])]>=4780/24/30, "survive","die")
  survive.mean.Pd[i] <- ifelse(df.bat1$Mean.Pd[which(df.bat1$Temp==df.data$T.a.round[i] & df.bat1$RH==df.data$RH.round[i])]>=4780/24/30, "survive","die")
  survive.high.Pd[i] <- ifelse(df.bat1$CIHigh.Pd[which(df.bat1$Temp==df.data$T.a.round[i] & df.bat1$RH==df.data$RH.round[i])]>=4780/24/30, "survive","die")
}
df.data <- cbind(df.data,survive.low,survive.mean,survive.high,survive.low.Pd,survive.mean.Pd,survive.high.Pd)


p1.l1 <- ggplot(df.bat1,aes(Temp,RH)) +
  theme_classic() + xlab("") + ylab("Relative Humidity (without Pd)") +
  geom_point(data=df.data, aes(T.a.round,RH.round,size=bout.length)) +
  scale_size_continuous(range=c(1,4),name="Measured torpor\nlength (days)") +
    theme(legend.justification=c(0,0.8))
p1.l2 <- ggplot(df.bat1,aes(Temp,RH)) +
  theme_classic() + xlab("") + ylab("Relative Humidity (without Pd)") +
  geom_raster(aes(fill=Mean.Pd)) + ggtitle("0.025 Quantile") +
  geom_point(data=df.data, aes(T.a.round,RH.round,shape=survive.mean.Pd)) +
  scale_shape_manual(values=c(7,16),name="Forecasted survival") +
  scale_fill_gradient2(low="#833437",mid="#8F8093",high="#67B9E9",
                       midpoint=5000/24/30, limits=c(0,12),
                       breaks=c(1,6,11),
                       labels=c("Mortality","Approximate winter length","Survival"),
                       name="Months to\nexhaustion") +
  theme(legend.justification=c(0,0.8))

l1 <- get_legend(p1.l1)
l2 <- get_legend(p1.l2)



p1 <- ggplot(df.bat1,aes(Temp,RH)) +
  theme_classic() + xlab("") + ylab("Relative Humidity (without Pd)") +
  geom_raster(aes(fill=CILow)) + ggtitle("0.025 Quantile") +
  geom_point(data=df.data, aes(T.a.round,RH.round,shape=survive.low, size=bout.length)) +
  scale_shape_manual(values=c(7,16),name="Forecasted survival") +
  scale_size_continuous(range=c(1,4),name="Measured torpor\nlength (days)") +
  scale_fill_gradient2(low="#833437",mid="#8F8093",high="#67B9E9",
                       midpoint=5000/24/30, limits=c(0,12),
                       breaks=c(1,6,11),
                       labels=c("Mortality","Approximate winter length","Survival"),
                       name="Months to\nexhaustion") +
  theme(legend.position="none")
p2 <- ggplot(df.bat1,aes(Temp,RH)) +
  theme_classic() + xlab("") + ylab("") +
  geom_raster(aes(fill=Mean)) + ggtitle("Mean") +
  geom_point(data=df.data, aes(T.a.round,RH.round,shape=survive.mean, size=bout.length)) +
  scale_shape_manual(values=c(7,16),name="Estimated survival") +
  scale_size_continuous(range=c(1,4),name="Measured\ntorpor length") +
  scale_fill_gradient2(low="#833437",mid="#8F8093",high="#67B9E9",midpoint=5000/24/30) +
  theme(legend.position="none")
p3 <- ggplot(df.bat1,aes(Temp,RH)) +
  theme_classic() + xlab("") + ylab("") +
  geom_raster(aes(fill=CIHigh)) + ggtitle("0.975 Quantile") +
  geom_point(data=df.data, aes(T.a.round,RH.round,shape=survive.high, size=bout.length)) +
  scale_shape_manual(values=c(16),name="Estimated survival") +
  scale_size_continuous(range=c(1,4),name="Measured\ntorpor length") +
  scale_fill_gradient2(low="#833437",mid="#8F8093",high="#67B9E9",midpoint=5000/24/30) +
  theme(legend.position="none")
p4 <- ggplot(df.bat1,aes(Temp,RH)) +
  theme_classic() + xlab("Temperature") + ylab("Relative Humidity (with Pd)") +
  geom_raster(aes(fill=CILow.Pd)) +
  geom_point(data=df.data, aes(T.a.round,RH.round,shape=survive.low.Pd, size=bout.length)) +
  scale_shape_manual(values=c(7,16),name="Estimated survival") +
  scale_size_continuous(range=c(1,4),name="Measured\ntorpor length") +
  scale_fill_gradient2(low="#833437",mid="#8F8093",high="#67B9E9",midpoint=5000/24/30) +
  theme(legend.position="none")
p5 <- ggplot(df.bat1,aes(Temp,RH)) +
  theme_classic() + xlab("Temperature") +  ylab("") +
  geom_raster(aes(fill=Mean.Pd)) +
  geom_point(data=df.data, aes(T.a.round,RH.round,shape=survive.mean.Pd, size=bout.length)) +
  scale_shape_manual(values=c(7,16),name="Estimated survival") +
  scale_size_continuous(range=c(1,4),name="Measured\ntorpor length") +
  scale_fill_gradient2(low="#833437",mid="#8F8093",high="#67B9E9",midpoint=5000/24/30) +
  theme(legend.position="none")
p6 <- ggplot(df.bat1,aes(Temp,RH)) +
  theme_classic() + xlab("Temperature") +  ylab("") +
  geom_raster(aes(fill=CIHigh.Pd)) +
  geom_point(data=df.data, aes(T.a.round,RH.round,shape=survive.high.Pd, size=bout.length)) +
  scale_shape_manual(values=c(7,16),name="Estimated survival") +
  scale_size_continuous(range=c(1,4),name="Measured\ntorpor length") +
  scale_fill_gradient2(low="#833437",mid="#8F8093",high="#67B9E9",midpoint=5000/24/30) +
  theme(legend.position="none")
# 
# legend <- get_legend(p1)
# p1 <- p1+theme(legend.position="none")
plot_grid(p1,p2,p3,l2,p4,p5,p6,l1,ncol=4)
# plot_grid(p1,p2,p3,legend,ncol=4)

p.survive.mean <- length(which(df.data$survive.mean=="survive"))/length(df.data$survive.mean)
p.survive.mean.Pd <- length(which(df.data$survive.mean.Pd=="survive"))/length(df.data$survive.mean.Pd)
p.survive.high <- length(which(df.data$survive.high=="survive"))/length(df.data$survive.high)
p.survive.high.Pd <- length(which(df.data$survive.high.Pd=="survive"))/length(df.data$survive.high.Pd)
