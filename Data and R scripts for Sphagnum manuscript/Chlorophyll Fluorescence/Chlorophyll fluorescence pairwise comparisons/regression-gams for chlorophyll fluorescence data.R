#gams for fluorescence data 
#chris goodall 02/08/22
#updated with linear models 02/20/22
#round by 50 -- 06-03-2023
#check the optimization -- 07-19-2023

#packages 
library(readr)
library(mgcv)
library(MASS)
library(reshape)

#import data
d1 <- read_csv("data/moss chl flu index by lat 04-14-2023.csv")

#clean data
# colnames(d1)
d1$lat <- round(d1$lat, 0) #round by lat
colnames(d1)[2] <- "Latitude" #change name 
d1$PPFD <- round_any(d1$PPFD, 50) #round to nearest increment of 50 for PPFD

#Subsetting to 725 PAR 
d1 <- d1[d1$PPFD > 749 & d1$PPFD < 751,] #subsetting to 750. 

#cat(paste(colnames(d1), collapse=', '))

#function to extract linear regression output values
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  s <- summary(modelobject)$coefficients[2,1] #slope
  f <- summary(modelobject)$fstatistic #f
  r <- summary(modelobject)$adj.r.squared #r^2 
  p <- pf(f[1],f[2],f[3],lower.tail=F) #p 
  attributes(p) <- NULL
  attributes(f) <- NULL
  f <- summary(modelobject)$fstatistic[1] #f again, but this time it's cool
  out <- cbind(s, r, f, p)
  return(out)
}


#linear models
FvFm_mean_lm <- lm(FvFm_mean ~ Latitude, data=d1)
a <- lmp(FvFm_mean_lm)
qP_lm <- lm(qP ~ Latitude, data=d1)
b <- lmp(qP_lm)
qL_lm <- lm(qL ~ Latitude, data=d1)
c<- lmp(qL_lm)
NPQ_lm <- lm(NPQ ~ Latitude, data=d1)
d<- lmp(NPQ_lm)
phiPSII_lm <- lm(phiPSII ~ Latitude, data=d1)
e <- lmp(phiPSII_lm)
phiNPQ_lm <- lm(phiNPQ ~ Latitude, data=d1)
f <- lmp(phiNPQ_lm)
phiNO_lm <- lm(phiNO ~ Latitude, data=d1)
g<- lmp(phiNO_lm)

#compile
(lm <- rbind(a,b,c,d,e,f,g))

#design function to extract values from GAM
lmp.2 <- function (modelobject) {
  dev <- summary(modelobject)[14] #explained deviance
  r2 <- summary(modelobject)$r.sq #r.squared
  edf <- summary(modelobject)$edf #EDF
  aic <- AIC(modelobject) #AIC
  p<- summary(modelobject)$s.table[,4] #p
  out <- cbind(dev,r2,edf,aic,p)
  return(out)
}

#use these to check gam
par(mfrow = c(2, 2))
#set seed for reproducibility
set.seed(26)

#make GAMs, check, store output
#FvFm
FvFm_mean_gam <- gam(FvFm_mean ~ s(Latitude, k=9), select = T, method = "REML", data=d1)
summary(FvFm_mean_gam)
gam.check(FvFm_mean_gam, rep=500) 
#plot(FvFm_mean_gam, residuals=TRUE,pch=19,cex=.3)
(a <- lmp.2(FvFm_mean_gam)) #fits.

#qP
qP_gam <- gam(qP ~ s(log(Latitude), k=9), method = "REML", data=d1)
summary(qP_gam)
gam.check(qP_gam, rep=500) #there seemed to be a log shaped curve spilling into the residuals, so I took the log of latitude to fit the model. it seemed to help?
#plot(qP_gam, residuals=TRUE,pch=19,cex=.3)
(b<- lmp.2(qP_gam)) #fits. 

#qL
qL_gam <- gam(qL ~ s(log(Latitude), k=9), method = "REML", data=d1)
summary(qL_gam)
gam.check(qL_gam) #similar to qP, I took the log of latitude here and it seems to have helped the model. 
#plot(qL_gam, residuals=TRUE,pch=19,cex=.3)
(c <- lmp.2(qL_gam)) #fits.

#NPQ
NPQ_gam <- gam(NPQ ~ s(Latitude, k=9), method = "REML", data=d1)
summary(NPQ_gam)
gam.check(NPQ_gam, rep=500)
#plot(NPQ_gam, residuals=TRUE,pch=19,cex=.3)
(d <- lmp.2(NPQ_gam)) #fits.

#phiPSII
phiPSII_gam <- gam(phiPSII ~ s(Latitude, k=9), select =T, method = "REML", data=d1)
summary(phiPSII_gam)
gam.check(phiPSII_gam, rep=500) 
#plot(phiPSII_gam, residuals = T, pch=19, cex=0.3) #fits.
(e <- lmp.2(phiPSII_gam))

#phiNPQ
phiNPQ_gam <- gam(phiNPQ ~ s(Latitude, k=9), select = T, method = "REML", data=d1)
summary(phiNPQ_gam)
gam.check(phiNPQ_gam, rep=500) #passes
#plot(phiNPQ_gam, residuals = T, pch=19, cex=0.3) #fits.
(f <- lmp.2(phiNPQ_gam))

#phiNO
phiNO_gam <- gam(phiNO_new ~ s(Latitude, k=9), method = "REML", data=d1)
summary(phiNO_gam)
gam.check(phiNO_gam, rep=500) #passes
#plot(phiNO_gam, residuals = T, pch=19, cex=0.3) #fits.
(g <- lmp.2(phiNO_gam))

#combine
(gam <- rbind(a,b,c,d,e,f,g))

#combine linear regression and GAM results
ID <- c("FvFm","qP","qL","NPQ","phiPSII","phiNPQ","phiNO")
(tab.1 <- cbind(ID, gam))
tab <- cbind(ID, lm, gam)

#rename columns 
colnames(tab) <- c("ID","slope","adj.r.squared","f","p","explained deviance","r.squared","EDF","AIC","p.gam")
tab

#write output
#write.csv(tab, "lmgam_chl.flu 04-15-2023.csv")
dev.off()


##

