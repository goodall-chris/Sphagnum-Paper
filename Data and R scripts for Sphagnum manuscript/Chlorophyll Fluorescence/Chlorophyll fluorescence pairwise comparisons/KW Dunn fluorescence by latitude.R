#KW/Dunn on Chl fluorescence -- latitude -- quantum yield
#chris goodall -- 02/02/22
#updated to have FvFm (dark adapted) 04-14-2023
#updated to change taxa factor order 04-23-2023
#round by 50 -- 06-03-2023
#cleaned code more 08-03-2023

#packages 
library(stats)
library(multcomp) #for analysis 
library(viridis) #for plotting
library(ggplot2)
library(dplyr)
library(stringr)
library(rstatix)
library(ggpubr)
library(tidyverse)
library(multcompView)
library(reshape)
library(readr)

#import data
d1 <- read_csv("data/chlorophyll fluorescence by latitude 08-31-2023.csv")

#clean data
d1$Latitude <- round(d1$Latitude, 0) #round by latitude
d1$PPFD <- round_any(d1$PPFD, 50) #round to nearest increment of 50 for PPFD
#d1$Latitude=as.factor(d1$Latitude)

#Subsetting to 750 PAR 
d1 <- d1[d1$PPFD > 749 & d1$PPFD < 751,] #subsetting data to 750 

#cat(paste0('c("', paste(colnames(d1), collapse='", "'), '")'))

#run nonparametric test: FvFm_mean
(res.kruskal <- d1 %>% kruskal_test(FvFm_mean ~ Latitude)) #this will indicate if there are significant differences between groups
d1 %>% kruskal_effsize(FvFm_mean ~ Latitude) #this indicates the amount of variance between groups explained by the treatment. 
#The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
#value of .482 is very large it suggests that latitudinal sites vary by FvFm

#Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
pwc <- d1 %>% 
    dunn_test(FvFm_mean ~ Latitude, p.adjust.method = "bonferroni") #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 

#designate significant differences between groups with letters
holder <- pwc$p #extracts p.values from dunn_test
names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
hold <-multcompLetters(holder) #this generates the letters 
letters.FvFm_mean <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector
letters.FvFm_mean

#create data frame for putting on letters
value_max_FvFm_mean = d1 %>%
  group_by(Latitude) %>%
  summarize(max_value = max(FvFm_mean))
value_max_FvFm_mean

# value_max_FvFm_mean = tapply(d1$FvFm_mean, d1$Latitude, max) %>% as.data.frame(.)
# value_max_FvFm_mean$Latitude <- rownames(value_max_FvFm_mean) %>% as.numeric(.)
# colnames(value_max_FvFm_mean)[1]="max_value"

#create text frame for graphic overlay
text_FvFm_mean <- cbind(value_max_FvFm_mean, letters.FvFm_mean) # combine both types of data 
text_FvFm_mean=text_FvFm_mean[,c(1,3,2)] #rearrange columns

#add sneaky features for plotting
d1$lat <- as.factor(d1$Latitude)
text_FvFm_mean$lat <- as.factor(text_FvFm_mean$Latitude) 

#color palette
pal.2 <- viridis(9, direction = 1) #create palette

#make plot
FvFm_mean <- ggplot(d1, aes(x = Latitude, y = FvFm_mean,
                            fill = lat)) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65, position="jitter")+
  geom_violin(alpha = 0.3, width =1) +
  scale_y_continuous(limits = c(min(d1$FvFm_mean)-.05, max(d1$FvFm_mean)+.05)) +
  theme_bw() +
  theme(text = element_text(size = 24),
        legend.position = "none",
        legend.title = element_blank()) +
  scale_fill_manual(values = pal.2) +
  geom_text(data = text_FvFm_mean, 
            aes(x=Latitude, y = 0.02 + max_value, label = value), 
            vjust=-0,
            size=4.5, #CHANGE AND NUDGE
            fontface = "bold") +
  labs(x= "Latitude", y=expression("F"[v]*"/F"[m]*""))
FvFm_mean

#run nonparametric test: qP
(res.kruskal <- d1 %>% kruskal_test(qP ~ Latitude)) #this will indicate if there are siginficant differences between groups
d1 %>% kruskal_effsize(qP ~ Latitude) #this indicates the amount of variance between groups explained by the treatment. 
#The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
#value of .431 is very large it suggests that latitudinal sites vary by qP 

#Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
(pwc <- d1 %>% 
    dunn_test(qP ~ Latitude, p.adjust.method = "bonferroni")) #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 

#designate significant differences between groups with letters
holder <- pwc$p #extracts p.values from dunn_test
names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
hold <-multcompLetters(holder) #this generates the letters 
letters.qP <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector
letters.qP

#create data frame for putting on letters
value_max_qP = d1 %>%
  group_by(Latitude) %>%
  summarize(max_value = max(qP))
value_max_qP

#create text frame for graphic overlay
text_qP <- cbind(value_max_qP, letters.qP) # combine both types of data 
text_qP=text_qP[,c(1,3,2)] #rearrange columns

#add sneaky features for plotting
text_qP$lat <- as.factor(text_qP$Latitude) 

#make plot
qP <- ggplot(d1, aes(x = Latitude, y = qP,
                     fill = lat)) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65, position="jitter")+
  geom_violin(alpha = 0.3, width =1) +
  scale_y_continuous(limits = c(min(d1$qP)-.05, max(d1$qP)+.05)) +
  theme_bw() +
  theme(text = element_text(size = 24),
        legend.position = "none",
        legend.title = element_blank()) +
  scale_fill_manual(values = pal.2) +
  geom_text(data = text_qP, 
            aes(x=Latitude, y = 0.02 + max_value, label = value), 
            vjust=.5,
            size=4.5, 
            fontface = "bold") 
qP


#run nonparametric test: qL
(res.kruskal <- d1 %>% kruskal_test(qL ~ Latitude)) #this will indicate if there are siginficant differences between groups
d1 %>% kruskal_effsize(qL ~ Latitude) #this indicates the amount of variance between groups explained by the treatment. 
#The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
#value of .376 is very large it suggests that latitudinal sites vary by qL 

#Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
(pwc <- d1 %>% 
    dunn_test(qL ~ Latitude, p.adjust.method = "bonferroni")) #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 

#designate significant differences between groups with letters
holder <- pwc$p #extracts p.values from dunn_test
names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
hold <-multcompLetters(holder) #this generates the letters 
letters.qL <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector
letters.qL

#create data frame for putting on letters
value_max_qL = d1 %>%
  group_by(Latitude) %>%
  summarize(max_value = max(qL))
value_max_qL

#create text frame for graphic overlay
text_qL <- cbind(value_max_qL, letters.qL) # combine both types of data 
text_qL=text_qL[,c(1,3,2)] #rearrange columns

#add sneaky features for plotting
text_qL$lat <- as.factor(text_qL$Latitude) 

#make plot
qL <- ggplot(d1, aes(x = Latitude, y = qL,
                     fill = lat)) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65,position="jitter")+
  geom_violin(alpha = 0.3, width =1) +
  scale_y_continuous(limits = c(min(d1$qL)-.05, max(d1$qL)+.05)) +
  theme_bw() +
  theme(text = element_text(size = 24),
        legend.position = "none",
        legend.title = element_blank()) +
  scale_fill_manual(values = pal.2) +
  geom_text(data = text_qL, 
            aes(x=Latitude, y = 0.02 + max_value, label = value), 
            vjust=-0.3,
            size=4.5, 
            fontface = "bold") 
qL


#run nonparametric test: NPQ
(res.kruskal <- d1 %>% kruskal_test(NPQ ~ Latitude)) #this will indicate if there are siginficant differences between groups
d1 %>% kruskal_effsize(NPQ ~ Latitude) #this indicates the amount of variance between groups explained by the treatment. 
#The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
#value of .172 is very large it suggests that latitudinal sites vary by NPQ 

#Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
(pwc <- d1 %>% 
    dunn_test(NPQ ~ Latitude, p.adjust.method = "bonferroni")) #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 

#designate significant differences between groups with letters
holder <- pwc$p #extracts p.values from dunn_test
names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
hold <-multcompLetters(holder) #this generates the letters 
letters.NPQ <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector
letters.NPQ

#create data frame for putting on letters
value_max_NPQ = d1 %>%
  group_by(Latitude) %>%
  summarize(max_value = max(NPQ))
value_max_NPQ

#create text frame for graphic overlay
text_NPQ <- cbind(value_max_NPQ, letters.NPQ) # combine both types of data 
text_NPQ=text_NPQ[,c(1,3,2)] #rearrange columns

#add sneaky features for plotting
text_NPQ$lat <- as.factor(text_NPQ$Latitude) 

#make plot
NPQ <- ggplot(d1, aes(x = Latitude, y = NPQ,
                      fill = lat)) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65, position="jitter")+
  geom_violin(alpha = 0.3, width =1) +
  scale_y_continuous(limits = c(min(d1$NPQ)-.1, max(d1$NPQ)+.15)) +
  theme_bw() +
  theme(text = element_text(size = 24),
        legend.position = "none",
        legend.title = element_blank()) +
  scale_fill_manual(values = pal.2) +
  geom_text(data = text_NPQ, 
            aes(x=Latitude, y = 0.05 + max_value, label = value), 
            vjust=-0.3,
            size=4.5, 
            fontface = "bold") 
NPQ


#run nonparametric test: phiPSII
(res.kruskal <- d1 %>% kruskal_test(phiPSII ~ Latitude)) #this will indicate if there are siginficant differences between groups
d1 %>% kruskal_effsize(phiPSII ~ Latitude) #this indicates the amount of variance between groups explained by the treatment. 
#The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
#value of .586 is very large it suggests that latitudinal sites vary by phiPSII 

#Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
(pwc <- d1 %>% 
    dunn_test(phiPSII ~ Latitude, p.adjust.method = "bonferroni")) #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 

#designate significant differences between groups with letters
holder <- pwc$p #extracts p.values from dunn_test
names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
hold <-multcompLetters(holder) #this generates the letters 
letters.phiPSII <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector
letters.phiPSII

#create data frame for putting on letters
value_max_phiPSII = d1 %>%
  group_by(Latitude) %>%
  summarize(max_value = max(phiPSII))
value_max_phiPSII

#create text frame for graphic overlay
text_phiPSII <- cbind(value_max_phiPSII, letters.phiPSII) # combine both types of data 
text_phiPSII=text_phiPSII[,c(1,3,2)] #rearrange columns

#add sneaky features for plotting
text_phiPSII$lat <- as.factor(text_phiPSII$Latitude) 

#color palette
pal.2 <- viridis(9, direction = 1) #create palette

#make plot
phiPSII <- ggplot(d1, aes(x = Latitude, y = phiPSII,
                          fill = lat)) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65,  position="jitter")+
  geom_violin(alpha = 0.3, width =1) +
  scale_y_continuous(limits = c(min(d1$phiPSII)-.04, max(d1$phiPSII)+.04)) +
  theme_bw() +
  theme(text = element_text(size = 24),
        legend.position = "none",
        legend.title = element_blank()) +
  scale_fill_manual(values = pal.2) +
  geom_text(data = text_phiPSII, 
            aes(x=Latitude, y = 0.01 + max_value, label = value), 
            vjust=-0.3,
            size=4.5, 
            fontface = "bold") +
  labs(y=expression("Φ"[PSII]))
phiPSII


#run nonparametric test: phiNPQ
(res.kruskal <- d1 %>% kruskal_test(phiNPQ ~ Latitude)) #this will indicate if there are siginficant differences between groups
d1 %>% kruskal_effsize(phiNPQ ~ Latitude) #this indicates the amount of variance between groups explained by the treatment. 
#The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
#value of .153 is very large it suggests that latitudinal sites vary by phiNPQ 

#Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
(pwc <- d1 %>% 
    dunn_test(phiNPQ ~ Latitude, p.adjust.method = "bonferroni")) #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 

#designate significant differences between groups with letters
holder <- pwc$p #extracts p.values from dunn_test
names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
hold <-multcompLetters(holder) #this generates the letters 
letters.phiNPQ <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector
letters.phiNPQ

#create data frame for putting on letters
value_max_phiNPQ = d1 %>%
  group_by(Latitude) %>%
  summarize(max_value = max(phiNPQ))
value_max_phiNPQ

#create text frame for graphic overlay
text_phiNPQ <- cbind(value_max_phiNPQ, letters.phiNPQ) # combine both types of data 
text_phiNPQ=text_phiNPQ[,c(1,3,2)] #rearrange columns

#add sneaky features for plotting
text_phiNPQ$lat <- as.factor(text_phiNPQ$Latitude) 

#make plot
phiNPQ <- ggplot(d1, aes(x = Latitude, y = phiNPQ,
                         fill = lat)) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65,  position="jitter")+
  geom_violin(alpha = 0.3, width =1) +
  scale_y_continuous(limits = c(min(d1$phiNPQ)-.05, max(d1$phiNPQ)+.05)) +
  theme_bw() +
  theme(text = element_text(size = 24),
        legend.position = "none",
        legend.title = element_blank()) +
  scale_fill_manual(values = pal.2) +
  geom_text(data = text_phiNPQ, 
            aes(x=Latitude, y = 0.01 + max_value, label = value), 
            vjust=-0.3,
            size=4.5, 
            fontface = "bold") +
  labs(y=expression("Φ"[NPQ]))
phiNPQ


#run nonparametric test: phiNO_new
(res.kruskal <- d1 %>% kruskal_test(phiNO_new ~ Latitude)) #this will indicate if there are siginficant differences between groups
d1 %>% kruskal_effsize(phiNO_new ~ Latitude) #this indicates the amount of variance between groups explained by the treatment. 
#The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
#value of .341 is very large it suggests that latitudinal sites vary by phiNO 

#Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
(pwc <- d1 %>% 
    dunn_test(phiNO_new ~ Latitude, p.adjust.method = "bonferroni")) #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 

#designate significant differences between groups with letters
holder <- pwc$p #extracts p.values from dunn_test
names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
hold <-multcompLetters(holder) #this generates the letters 
letters.phiNO_new <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector
letters.phiNO_new

#create data frame for putting on letters
value_max_phiNO_new = d1 %>%
  group_by(Latitude) %>%
  summarize(max_value = max(phiNO_new))
value_max_phiNO_new

#create text frame for graphic overlay
text_phiNO_new <- cbind(value_max_phiNO_new, letters.phiNO_new) # combine both types of data 
text_phiNO_new=text_phiNO_new[,c(1,3,2)] #rearrange columns

#add sneaky features for plotting
text_phiNO_new$lat <- as.factor(text_phiNO_new$Latitude) 

#make plot
phiNO <- ggplot(d1, aes(x = Latitude, y = phiNO_new,
                        fill = lat)) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65,  position="jitter")+
  geom_violin(alpha = 0.3, width =1) +
  scale_y_continuous(limits = c(min(d1$phiNO_new)-.05, max(d1$phiNO_new)+.05)) +
  theme_bw() +
  theme(text = element_text(size = 24),
        legend.position = "none",
        legend.title = element_blank()) +
  scale_fill_manual(values = pal.2) +
  geom_text(data = text_phiNO_new, 
            aes(x=Latitude, y = 0.01 + max_value, label = value), 
            vjust=-0.3,
            size=4.5, 
            fontface = "bold") +
  labs(y=expression("Φ"[NO]))
phiNO


#make fake plot for legend
library(cowplot)
library(grid)
library(gridExtra) 

d1$Latitude=as.factor(d1$Latitude)
filler <- ggplot(d1, aes(x = Latitude, y = FvFm_mean,
                         fill = Latitude)) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65) +
  geom_violin(alpha = 0.3, width =2) +
  scale_y_continuous(limits = c(min(d1$FvFm_mean)-.05, max(d1$FvFm_mean)+.05)) +
  theme_bw() +
  scale_fill_manual(values = pal.2) +
  theme(text = element_text(size = 24),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        legend.text=element_text(size=22),
        legend.title.align=0.5,
        legend.title=element_text(size=32, face="bold"),
        legend.position = "right") +
  guides(color=guide_legend(ncol=3))
legend <- cowplot::get_legend(filler) #pull legend 


#create layout and grid arrange plots 
library(patchwork)
layout <- '
ABCD
EFGX'

# tiff("chl flu -- latitude.tiff", units="in", width=20, height=10, res=300)
# wrap_plots(A = FvFm_mean, B = qP, C = qL, D = NPQ, E= phiPSII,
#            `F` = phiNPQ, G = phiNO, X = legend, design=layout) +
#   plot_annotation(tag_levels = 'A')
# dev.off()

