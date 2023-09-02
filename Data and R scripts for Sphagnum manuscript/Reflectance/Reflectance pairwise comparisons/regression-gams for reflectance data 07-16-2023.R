#trying GAMs
#chris goodall -- 02/06/22
#checking gam script -- 07-19-2023

#packages 
library(readr)
library(mgcv)
library(MASS)

#import data
d1 <- read.csv("data/reflectance by latitude.csv")

#clean data
#d1$Latitude <- round(d1$Latitude, 0) #round by latitutde

#function to extract linear regression output values
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  slope <- summary(modelobject)$coefficients[2,1] #slope
  r_squared <- summary(modelobject)$adj.r.squared #r^2 
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F) #p 
  attributes(p) <- NULL
  attributes(f) <- NULL
  f <- summary(modelobject)$fstatistic[1] #f again, but this time it's cool
  out <- cbind(slope, r_squared, f, p)
  return(out)
}

#Run linear models
NDVI_lm <- lm(NDVI ~ Latitude, data=d1)
(a<- lmp(NDVI_lm))
PRI_lm <- lm(PRI ~ Latitude, data=d1)
(b <- lmp(PRI_lm))
ARI_lm <- lm(ARI ~ Latitude, data=d1)
c <- lmp(ARI_lm)
ExGm_lm <- lm(ExGm ~ Latitude, data=d1)
d <- lmp(ExGm_lm)
GRVI_lm <- lm(GRVI ~ Latitude, data=d1)
e <- lmp(GRVI_lm)
RGR_lm <- lm(RGR ~ Latitude, data=d1)
f <- lmp(RGR_lm)

#compile lm outputs
(lm <- rbind(a,b,c,d,e,f))

#design function to extract values from GAM
lmp.2 <- function (modelobject) {
  dev <- summary(modelobject)[14] #explained deviance
  r2 <- summary(modelobject)$r.sq #r.squared
  edf <- summary(modelobject)$edf #EDF
  aic <- AIC(modelobject) #AIC
  p <- summary(modelobject)$s.table[,4] #p
  out <- cbind(dev,r2,edf,aic,p)
  return(out)
}

#check to see gam plots
par(mfrow=c(2,2))

#set seed for reproducibility
set.seed(26)

#make GAMs, check, store output
#NDVI
NDVI_gam <- gam(NDVI ~ s(Latitude, k=9), method = "REML", select = T, data=d1)
summary(NDVI_gam)
gam.check(NDVI_gam, rep=500) #passes
#plot(NDVI_gam, residuals=TRUE,pch=19,cex=.3)
(a<- lmp.2(NDVI_gam))


#and altitude
# NDVI_gam <- gam(NDVI ~ te(Altitude, Latitude), method = "REML", select = T, data=d1) #need to determine if tensor product makes sense.
# summary(NDVI_gam)
# gam.check(NDVI_gam, rep=500) #passes
# #plot(NDVI_gam, residuals=TRUE,pch=19,cex=.3)
##


#PRI
PRI_gam <- gam(PRI ~ s(Latitude, k=9), select = T, method = "REML", data=d1) 
summary(PRI_gam)
gam.check(PRI_gam, rep=500) #passes, there is a gap in the residuals vs. predictor plot though.
#plot(PRI_gam, residuals=TRUE,pch=19,cex=.3)
b <- lmp.2(PRI_gam)

#and altitude
# PRI_gam2 <- gam(PRI ~ te(Altitude, Latitude), select = T, method = "REML", data=d1)
# summary(PRI_gam2)
# gam.check(PRI_gam2, rep=500) #passes, there is a gap in the residuals vs. predictor plot though.
#plot(PRI_gam, residuals=TRUE,pch=19,cex=.3)

#ARI
ARI_gam <- gam(ARI ~ s(Latitude, k=9), method = "REML", data=d1) 
summary(ARI_gam)
gam.check(ARI_gam, rep=500) #passes.
#plot(ARI_gam, residuals=TRUE,pch=19,cex=.3)
c <- lmp.2(ARI_gam) 

#and altitude
# ARI_gam2 <- gam(ARI ~ te(Altitude, Latitude), method = "REML", data=d1)
# summary(ARI_gam2)
# gam.check(ARI_gam2, rep=500) #passes.
# #plot(ARI_gam, residuals=TRUE,pch=19,cex=.3)

#ExGm
ExGm_gam <- gam(ExGm ~ s(Latitude, k=9), select=T, method = "REML", data=d1)
summary(ExGm_gam)
gam.check(ExGm_gam, rep=500) #passes
#plot(ExGm_gam, residuals=TRUE,pch=19,cex=.3)
d <- lmp.2(ExGm_gam)

#and altitude
ExGm_gam2 <- gam(ExGm ~te(Altitude, Latitude), select=T, method = "REML", data=d1)
summary(ExGm_gam2)
gam.check(ExGm_gam, rep=500) #passes
#plot(ExGm_gam, residuals=TRUE,pch=19,cex=.3)

#GRVI
GRVI_gam <- gam(GRVI ~ s(Latitude, k=9), select=T, method = "REML", data=d1)
summary(GRVI_gam)
gam.check(GRVI_gam, rep=500) #passes
#plot(GRVI_gam, residuals=TRUE,pch=19,cex=.3)
e <- lmp.2(GRVI_gam)

#add altitude
# GRVI_gam2 <- gam(GRVI ~ te(Altitude, Latitude), select=T, method = "REML", data=d1)
# summary(GRVI_gam2)
# gam.check(GRVI_gam2, rep=500) #passes
# #plot(GRVI_gam, residuals=TRUE,pch=19,cex=.3)
# e <- lmp.2(GRVI_gam)

#RGR
RGR_gam <- gam(RGR ~ s(Latitude, k=9), select=T, method = "REML", data=d1)
summary(RGR_gam)
gam.check(RGR_gam, rep=500) #passes
#plot(RGR_gam, residuals=TRUE,pch=19,cex=.3)
f <- lmp.2(RGR_gam)

#and altitude
# RGR_gam2 <- gam(RGR ~ te(Altitude, Latitude), select=T, method = "REML", data=d1)
# summary(RGR_gam2)
# gam.check(RGR_gam2, rep=500) #passes
# #plot(RGR_gam, residuals=TRUE,pch=19,cex=.3)
# f <- lmp.2(RGR_gam)

#compile gam data
(gam <- rbind(a,b,c,d,e,f))

#combine

#combine linear regression and GAM results
ID <- c("NDVI", "PRI", "ARI", "ExGm", "GRVI", "RGR")
tab <- cbind(ID, lm, gam)

#rename columns 
colnames(tab) <- c("ID","slope","f","r.squared","p","explained deviance","r.squared","EDF","AIC","p")
tab

#write output
#write.csv(tab, "lmgam_ref.csv")

#cat(paste(colnames(d1), collapse='," '))
#dev.off()
##
