#Plotting fluorescence data by latitude
#chris goodall -- 09.14.21
#re-wrote -- 01/24/22 
#updated to include varying PAR plots -- 1/30/22
#first draft complete -- 02/02/22
#round by 50 -- 06-03-2023

#packages 
library(stats)
library(multcomp) #for analysis 
library(viridis) #for plotting
library(ggplot2)
library(dplyr)
library(agricolae)
library(stringr)
library(rstatix)
library(ggpubr)
library(tidyverse)
library(multcompView)
library(reshape)
library(readr)
library(plyr)

#import data
d1 <- read_csv("data/moss chl flu index by lat 04-14-2023.csv")

#clean data
d1$lat <- round(d1$lat, 0) #round by lat
colnames(d1)[2] <- "Latitude" #change name 
d1$PPFD <- round_any(d1$PPFD, 50) #round to nearest increment of 50 for PPFD
d1$PPFD= as.factor(d1$PPFD) #make PPFD as.factor for plotting

rev(unique(d1$PPFD))

#make plots for general trends 
#fv/fm
fvfm <- ggplot(d1, aes(x=Latitude, y=FvFm_mean)) +
  geom_point(shape=16, alpha = 0.6, aes(color=PPFD)) + 
  geom_smooth(aes(color=PPFD), size = .7, se=F) +
  theme_bw() +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(x=expression('Latitude'),
       y=expression("Fv/Fm")) +
  theme(text = element_text(size = 24),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18, face="bold"), legend.position = "none") + 
  labs(y=expression("F"[v]*"/F"[m]*""))
#qP
qP <- ggplot(d1, aes(x=Latitude, y=qP)) +
  geom_point(shape=16, alpha = 0.6, aes(color=PPFD)) + 
  geom_smooth(aes(color=PPFD), size = .7, se=F) +
  theme_bw() +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(x=expression('Latitude'),
       y=expression("qP")) +
  theme(text = element_text(size = 24),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18, face="bold"), legend.position = "none")
#qL
qL <- ggplot(d1, aes(x=Latitude, y=qL)) +
  geom_point(shape=16, alpha = 0.6, aes(color=PPFD)) + 
  geom_smooth(aes(color=PPFD), size = .7, se=F) +
  theme_bw() +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(x=expression('Latitude'),
       y=expression("qL")) +
  theme(text = element_text(size = 24),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18, face="bold"), legend.position = "none")
#qN
qN <- ggplot(d1, aes(x=Latitude, y=qN)) +
  geom_point(shape=16, alpha = 0.6, aes(color=PPFD)) + 
  geom_smooth(aes(color=PPFD), size = .7, se=F) +
  theme_bw() +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(x=expression('Latitude'),
       y=expression("qN")) +
  theme(text = element_text(size = 24),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18, face="bold"), legend.position = "none")
#NPQ
NPQ <- ggplot(d1, aes(x=Latitude, y=NPQ)) +
  geom_point(shape=16, alpha = 0.6, aes(color=PPFD)) + 
  geom_smooth(aes(color=PPFD), size = .7, se=F) +
  theme_bw() +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(x=expression('Latitude'),
       y=expression("NPQ")) +
  theme(text = element_text(size = 24),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18, face="bold"), legend.position = "none")
#phiPSII
phiPSII <- ggplot(d1, aes(x=Latitude, y=phiPSII)) +
  geom_point(shape=16, alpha = 0.6, aes(color=PPFD)) + 
  geom_smooth(aes(color=PPFD), size = .7, se=F) +
  theme_bw() +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(x=expression('Latitude'),
      y=expression("Φ"[PSII])) +
  theme(text = element_text(size = 24),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18, face="bold"), legend.position = "none")
#phiNPQ
phiNPQ <- ggplot(d1, aes(x=Latitude, y=phiNPQ)) +
  geom_point(shape=16, alpha = 0.6, aes(color=PPFD)) + 
  geom_smooth(aes(color=PPFD), size = .7, se=F) +
  theme_bw() +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(x=expression('Latitude'),
       y=expression("Φ"[NPQ])) +
  theme(text = element_text(size = 24),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18, face="bold"), legend.position = "none")
#phiNO_new
phiNO <- ggplot(d1, aes(x=Latitude, y=phiNO_new)) +
  geom_point(shape=16, alpha = 0.6, aes(color=PPFD)) + 
  geom_smooth(aes(color=PPFD), size = .7, se=F) +
  theme_bw() +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(x=expression('Latitude'),
       y=expression("Φ"[NO])) +
  theme(text = element_text(size = 24),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18, face="bold"), legend.position = "none")
#JETR
JETR <- ggplot(d1, aes(x=Latitude, y=JETR)) +
  geom_point(shape=16, alpha = 0.6, aes(color=PPFD)) + 
  geom_smooth(aes(color=PPFD), size = .7, se=F) +
  theme_bw() +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(x=expression('Latitude'),
       y=expression("J"[ETR])) +
  theme(text = element_text(size = 24),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18, face="bold"), legend.position = "none")
#JNPQ
JNPQ <- ggplot(d1, aes(x=Latitude, y=JNPQ)) +
  geom_point(shape=16, alpha = 0.6, aes(color=PPFD)) + 
  geom_smooth(aes(color=PPFD), size = .7, se=F) +
  theme_bw() +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(x=expression('Latitude'),
       y=expression("J"[NPQ])) +
  theme(text = element_text(size = 24),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18, face="bold"), legend.position = "none")
#JNO
JNO <- ggplot(d1, aes(x=Latitude, y=JNO)) +
  geom_point(shape=16, alpha = 0.6, aes(color=PPFD)) + 
  geom_smooth(aes(color=PPFD), size = .7, se=F) +
  theme_bw() +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(x=expression('Latitude'),
       y=expression("J"[NO])) +
  theme(text = element_text(size = 24),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18, face="bold"), legend.position = "none")
#JTOT
Jtot <- ggplot(d1, aes(x=Latitude, y=Jtot)) +
  geom_point(shape=16, alpha = 0.6, aes(color=PPFD)) + 
  geom_smooth(aes(color=PPFD), size = .7, se=F) +
  theme_bw() +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(x=expression('Latitude'),
       y=expression("J"[tot])) +
  theme(text = element_text(size = 24),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18, face="bold"), legend.position = "none")

#make fake plot for legend
library(cowplot)
library(grid)
library(gridExtra) 
filler <- ggplot(d1, aes(x=Latitude, y=FvFm_mean)) +
  geom_point(shape=16, alpha = 0.6, aes(color=PPFD)) + 
  geom_smooth(aes(color=PPFD), size = .7, se=F) +
  theme_bw() +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(x=expression('Latitude'),
       y=expression("Fv'/Fm'")) +
  theme(text = element_text(size = 24),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        legend.text=element_text(size=22),
        legend.title.align=0.5,
        legend.title=element_text(size=32, face="bold"),
        legend.position = "right") +
  guides(color=guide_legend(ncol=3))
filler
legend <- cowplot::get_legend(filler) #pull legend 


#create layout and grid arrange plots 
library(patchwork)
layout <- '
ABCD
EFGX
HIJK'

# tiff("test.tiff", units="in", width=22, height=15, res=300)
# wrap_plots(A = fvfm, B = qP, C = qL, D= NPQ, E = phiPSII, F = phiNPQ,
#            G = phiNO, X = legend, H = JETR, I =JNPQ, J =JNO,
#            K = Jtot, design=layout) + plot_annotation(tag_levels = 'A')
# dev.off()
# 
# ##


#create layout and grid arrange plots for just quantum yield
library(patchwork)
layout <- '
ABCD
EFGX'

tiff("chl flu latitude 06-03-2023.tiff", units="in", width=20, height=11.5, res=300)
wrap_plots(A = fvfm, B = qP, C = qL, D = NPQ, E = phiPSII,
           `F` = phiNPQ, G = phiNO, X = legend, design=layout) + 
           plot_annotation(tag_levels = 'A')
dev.off()

##
