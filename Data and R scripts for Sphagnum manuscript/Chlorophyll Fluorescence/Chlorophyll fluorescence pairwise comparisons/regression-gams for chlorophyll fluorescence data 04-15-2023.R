#gams for fluorescence data 
#chris goodall 02/08/22
#updated with linear models 02/20/22
#round by 50 -- 06-03-2023

#packages 
library(readr)
library(mgcv)
library(visibly)
library(reshape)

#import data
d1 <- read_csv("data/moss chl flu index by lat 04-14-2023.csv")
colnames(d1)
#clean data
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
qN_lm <- lm(qN ~ Latitude, data=d1)
d <- lmp(qN_lm)
NPQ_lm <- lm(NPQ ~ Latitude, data=d1)
e<- lmp(NPQ_lm)
phiPSII_lm <- lm(phiPSII ~ Latitude, data=d1)
f <- lmp(phiPSII_lm)
phiNPQ_lm <- lm(phiNPQ ~ Latitude, data=d1)
g <- lmp(phiNPQ_lm)
phiNO_lm <- lm(phiNO ~ Latitude, data=d1)
h<- lmp(phiNO_lm)
JETR_lm <- lm(JETR ~ Latitude, data=d1)
i <- lmp(JETR_lm)
JNPQ_lm <- lm(JNPQ ~ Latitude, data=d1)
j <- lmp(JNPQ_lm)
JNO_lm <- lm(JNO ~ Latitude, data=d1)
k <- lmp(JNO_lm)
Jtot_lm <- lm(Jtot ~ Latitude, data=d1)
l <- lmp(Jtot_lm)

#compile
(lm <- rbind(a,b,c,d,e,f,g,h,i,j,k,l))

#design function to extract values from GAM
lmp.2 <- function (modelobject) {
  d<- summary(modelobject)[14] #explained deviance
  r <- summary(modelobject)$r.sq #r.squared
  e<- summary(modelobject)$edf #EDF
  a<- AIC(modelobject) #AIC
  p<- summary(modelobject)$s.table[,4] #p
  out <- cbind(d,r,e,a,p)
  return(out)
}

#use these to check gam
par(mfrow = c(2, 2))
#set seed for reproducibility
set.seed(26)

#make GAMs, check, store output
FvFm_mean_gam <- gam(FvFm_mean ~ s(Latitude, bs = "cs", k=9), method = "REML", data=d1)
gam.check(FvFm_mean_gam) #passes
(a <- lmp.2(FvFm_mean_gam))

qP_gam <- gam(qP ~ s(Latitude, bs = "cc", k=9), method = "REML", data=d1)
gam.check(qP_gam)
(b<- lmp.2(qP_gam)) #barely fit well 

qL_gam <- gam(qL ~ s(Latitude, bs = "cs", k=9), method = "REML", data=d1)
gam.check(qL_gam)
(c <- lmp.2(qL_gam)) #is not fit well 

qN_gam <- gam(qN ~ s(Latitude, bs = "cs", k=9), method = "REML", data=d1)
gam.check(qN_gam)
(d <- lmp.2(qN_gam))

NPQ_gam <- gam(NPQ ~ s(Latitude, bs = "cs", k=8), method = "REML", data=d1)
gam.check(NPQ_gam)
(e <- lmp.2(NPQ_gam)) #passable

phiPSII_gam <- gam(phiPSII ~ s(Latitude, k=9, bs = "cs" ), method = "REML", data=d1)
gam.check(phiPSII_gam) #passes
(f <- lmp.2(phiPSII_gam))
f
phiNPQ_gam <- gam(phiNPQ ~ s(Latitude, k=9, bs = "cs"), method = "REML", data=d1)
gam.check(phiNPQ_gam) #passes
(g <- lmp.2(phiNPQ_gam))

phiNO_gam <- gam(phiNO ~ s(Latitude, bs = "cs", k=9), method = "REML", data=d1)
gam.check(phiNO_gam) #passes
(h <- lmp.2(phiNO_gam))

JETR_gam <- gam(JETR ~ s(Latitude, bs = "cc", k=9), method = "REML", data=d1)
gam.check(JETR_gam) #passes
(i <- lmp.2(JETR_gam))

JNPQ_gam <- gam(JNPQ ~ s(Latitude, bs = "bs", k=9), method = "REML", data=d1)
gam.check(JNPQ_gam) #passes
(j <- lmp.2(JNPQ_gam))

JNO_gam <- gam(JNO ~ s(Latitude, bs = "cs", k=9), method = "REML", data=d1)
gam.check(JNO_gam) #passes, tight call
k <- lmp.2(JNO_gam)

Jtot_gam <- gam(Jtot ~ s(Latitude, bs = "cs", k=9), method = "REML", data=d1)
gam.check(Jtot_gam) #also broken. 
l <- lmp.2(Jtot_gam)

#combine
(gam <- rbind(a,b,c,d,e,f,g,h,i,j,k,l))

#combine linear regression and GAM results
ID <- c("FvFm","qP","qL","qN", "NPQ","phiPSII","phiNPQ","phiNO","JETR","JNPQ","JNO","Jtot")
tab <- cbind(ID, lm, gam)

#rename columns 
colnames(tab) <- c("ID","slope","adj.r.squared","f","p","explained deviance","r.squared","EDF","AIC","p.gam")
tab

#write output
#write.csv(tab, "lmgam_chl.flu 04-15-2023.csv")
dev.off()


##

