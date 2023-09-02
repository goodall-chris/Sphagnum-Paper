#KW/Dunn on Chl fluorescence -- taxa
#chris goodall -- 02/02/22
#updated by taxa -- 11/06/2022
#updated to have FvFm (dark adapted) 04-14-2023
#updated to change taxa factor order 04-23-2023
#round by 50 -- 06-03-2023

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
d1 <- read_csv("data/chlorophyll fluorescence by taxa.csv")


#clean data
d1$PPFD <- round_any(d1$PPFD, 50) #round to nearest increment of 50 for PPFD

corn <- d1[d1$Species=="Zea",]

marchantia <- d1[d1$Species=="Marchantia",]

#Subsetting to 750 PAR 
d1 <- d1[d1$PPFD > 749 & d1$PPFD < 751,] #subsetting to 750 PPFD

d2 <- d1[!d1$Species=="Zea",]
#cat(paste0('c("', paste(colnames(d1), collapse='", "'), '")'))

d1$Species <- factor(d1$Species, levels = c("Sphagnum", "Marchantia", "Zea")) #re-order data so sphagnum is read first

#run nonparametric test: FvFm_mean
(res.kruskal <- d1 %>% kruskal_test(FvFm_mean ~ Species)) #this will indicate if there are siginficant differences between groups
d1 %>% kruskal_effsize(FvFm_mean ~ Species) #this indicates the amount of variance between groups explained by the treatment. 
#The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
#value of .486 is very large it suggests that latitudinal sites vary by FvFm_mean 

#Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
(pwc <- d1 %>% 
    dunn_test(FvFm_mean ~ Species, p.adjust.method = "bonferroni")) #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 

#designate significant differences between groups with letters
holder <- pwc$p #extracts p.values from dunn_test
names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
hold <-multcompLetters(holder) #this generates the letters 
letters.FvFm_mean <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector
letters.FvFm_mean

#create data frame for putting on letters
value_max_FvFm_mean = d1 %>%
  group_by(Species) %>%
  summarize(max_value = max(FvFm_mean))
value_max_FvFm_mean

#create text frame for graphic overlay
text_FvFm_mean <- cbind(value_max_FvFm_mean, letters.FvFm_mean) # combine both types of data 
text_FvFm_mean=text_FvFm_mean[,c(1,3,2)] #rearrange columns

#create color palette
pal.2 <- plasma(3)[c(2,1,3)]

#make plot
FvFm_mean <- ggplot(d1, aes(x = Species, y = FvFm_mean,
                            fill = Species)) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65, position = "jitter") +
  geom_violin(alpha = 0.3, width =1.25) +
  scale_y_continuous(limits = c(min(d1$FvFm_mean)-.05, max(d1$FvFm_mean)+.05)) +
  theme_bw() +
  theme(text = element_text(size = 24), 
        axis.text.x = element_text(angle = 30, vjust = .6, hjust=.5),
        legend.position = "none", axis.title.x=element_blank(),
        legend.title = element_blank()) +
  scale_fill_manual(values = pal.2) +
  geom_text(data = text_FvFm_mean, 
            aes(x=Species, y = 0.02 + max_value, label = value), 
            vjust=-0.3,
            size=7, 
            fontface = "bold") +
  labs(y=expression("Fv/Fm"))
FvFm_mean

#run nonparametric test: qP
(res.kruskal <- d1 %>% kruskal_test(qP ~ Species)) #this will indicate if there are siginficant differences between groups
d1 %>% kruskal_effsize(qP ~ Species) #this indicates the amount of variance between groups explained by the treatment. 
#The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
#value of .454 is very large it suggests that latitudinal sites vary by qP 

#Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
(pwc <- d1 %>% 
    dunn_test(qP ~ Species, p.adjust.method = "bonferroni")) #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 

#designate significant differences between groups with letters
holder <- pwc$p #extracts p.values from dunn_test
names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
hold <-multcompLetters(holder) #this generates the letters 
letters.qP <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector
letters.qP

#create data frame for putting on letters
value_max_qP = d1 %>%
  group_by(Species) %>%
  summarize(max_value = max(qP))
value_max_qP

#create text frame for graphic overlay
text_qP <- cbind(value_max_qP, letters.qP) # combine both types of data 
text_qP=text_qP[,c(1,3,2)] #rearrange columns




#make plot
qP <- ggplot(d1, aes(x = Species, y = qP,
                     fill = Species)) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65, position = "jitter") +
  geom_violin(alpha = 0.3, width =1.25) +
  scale_y_continuous(limits = c(min(d1$qP)-.05, max(d1$qP)+.05)) +
  theme_bw() +
  theme(text = element_text(size = 24), 
        axis.text.x = element_text(angle = 30, vjust = .6, hjust=.5),
        legend.position = "none", axis.title.x=element_blank(),
        legend.title = element_blank()) +
  scale_fill_manual(values = pal.2) +
  geom_text(data = text_qP, 
            aes(x=Species, y = 0.02 + max_value, label = value), 
            vjust=-0.3,
            size=7, 
            fontface = "bold") 

#run nonparametric test: qL
(res.kruskal <- d1 %>% kruskal_test(qL ~ Species)) #this will indicate if there are siginficant differences between groups
d1 %>% kruskal_effsize(qL ~ Species) #this indicates the amount of variance between groups explained by the treatment. 
#The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
#value of .530 is very large it suggests that latitudinal sites vary by qL 

#Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
(pwc <- d1 %>% 
    dunn_test(qL ~ Species, p.adjust.method = "bonferroni")) #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 

#designate significant differences between groups with letters
holder <- pwc$p #extracts p.values from dunn_test
names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
hold <-multcompLetters(holder) #this generates the letters 
letters.qL <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector
letters.qL

#create data frame for putting on letters
value_max_qL = d1 %>%
  group_by(Species) %>%
  summarize(max_value = max(qL))
value_max_qL

#create text frame for graphic overlay
text_qL <- cbind(value_max_qL, letters.qL) # combine both types of data 
text_qL=text_qL[,c(1,3,2)] #rearrange columns




#make plot
qL <- ggplot(d1, aes(x = Species, y = qL,
                     fill = Species)) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65, position = "jitter") +
  geom_violin(alpha = 0.3, width =1.25) +
  scale_y_continuous(limits = c(min(d1$qL)-.05, max(d1$qL)+.05)) +
  theme_bw() +
  theme(text = element_text(size = 24), 
        axis.text.x = element_text(angle = 30, vjust = .6, hjust=.5),
        legend.position = "none", axis.title.x=element_blank(),
        legend.title = element_blank()) +
  scale_fill_manual(values = pal.2) +
  geom_text(data = text_qL, 
            aes(x=Species, y = 0.02 + max_value, label = value), 
            vjust=-0.3,
            size=7, 
            fontface = "bold") 

#run nonparametric test: NPQ
(res.kruskal <- d1 %>% kruskal_test(NPQ ~ Species)) #this will indicate if there are siginficant differences between groups
d1 %>% kruskal_effsize(NPQ ~ Species) #this indicates the amount of variance between groups explained by the treatment. 
#The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
#value of .535 is very large it suggests that latitudinal sites vary by NPQ 

#Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
(pwc <- d1 %>% 
    dunn_test(NPQ ~ Species, p.adjust.method = "bonferroni")) #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 

#designate significant differences between groups with letters
holder <- pwc$p #extracts p.values from dunn_test
names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
hold <-multcompLetters(holder) #this generates the letters 
letters.NPQ <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector
letters.NPQ

#create data frame for putting on letters
value_max_NPQ = d1 %>%
  group_by(Species) %>%
  summarize(max_value = max(NPQ))
value_max_NPQ

#create text frame for graphic overlay
text_NPQ <- cbind(value_max_NPQ, letters.NPQ) # combine both types of data 
text_NPQ=text_NPQ[,c(1,3,2)] #rearrange columns




#make plot
NPQ <- ggplot(d1, aes(x = Species, y = NPQ,
                      fill = Species)) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65, position = "jitter") +
  geom_violin(alpha = 0.3, width =1.25) +
  scale_y_continuous(limits = c(min(d1$NPQ)-.1, max(d1$NPQ)+.4)) +
  theme_bw() +
  theme(text = element_text(size = 24), axis.text.x = element_text(angle = 30, vjust = .6, hjust=.5),
        legend.position = "none", axis.title.x=element_blank(),
        legend.title = element_blank()) +
  scale_fill_manual(values = pal.2) +
  geom_text(data = text_NPQ, 
            aes(x=Species, y = 0.1 + max_value, label = value), 
            vjust=-0.3,
            size=7, 
            fontface = "bold") 

#run nonparametric test: phiPSII
(res.kruskal <- d1 %>% kruskal_test(phiPSII ~ Species)) #this will indicate if there are siginficant differences between groups
d1 %>% kruskal_effsize(phiPSII ~ Species) #this indicates the amount of variance between groups explained by the treatment. 
#The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
#value of .468 is very large it suggests that latitudinal sites vary by phiPSII 

#Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
(pwc <- d1 %>% 
    dunn_test(phiPSII ~ Species, p.adjust.method = "bonferroni")) #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 

#designate significant differences between groups with letters
holder <- pwc$p #extracts p.values from dunn_test
names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
hold <-multcompLetters(holder) #this generates the letters 
letters.phiPSII <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector
letters.phiPSII

#create data frame for putting on letters
value_max_phiPSII = d1 %>%
  group_by(Species) %>%
  summarize(max_value = max(phiPSII))
value_max_phiPSII

#create text frame for graphic overlay
text_phiPSII <- cbind(value_max_phiPSII, letters.phiPSII) # combine both types of data 
text_phiPSII=text_phiPSII[,c(1,3,2)] #rearrange columns




#make plot
phiPSII <- ggplot(d1, aes(x = Species, y = phiPSII,
                          fill = Species)) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65, position = "jitter") +
  geom_violin(alpha = 0.3, width =1.25) +
  scale_y_continuous(limits = c(min(d1$phiPSII)-.05, max(d1$phiPSII)+.05)) +
  theme_bw() +
  theme(text = element_text(size = 24), axis.text.x = element_text(angle = 30, vjust = .6, hjust=.5),
        legend.position = "none", axis.title.x=element_blank(),
        legend.title = element_blank()) +
  scale_fill_manual(values = pal.2) +
  geom_text(data = text_phiPSII, 
            aes(x=Species, y = 0.01 + max_value, label = value), 
            vjust=-0.3,
            size=7, 
            fontface = "bold") +
  labs(y=expression("Φ"[PSII]))


#run nonparametric test: phiNPQ
(res.kruskal <- d1 %>% kruskal_test(phiNPQ ~ Species)) #this will indicate if there are siginficant differences between groups
d1 %>% kruskal_effsize(phiNPQ ~ Species) #this indicates the amount of variance between groups explained by the treatment. 
#The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
#value of .540 is very large it suggests that latitudinal sites vary by phiNPQ 

#Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
(pwc <- d1 %>% 
    dunn_test(phiNPQ ~ Species, p.adjust.method = "bonferroni")) #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 

#designate significant differences between groups with letters
holder <- pwc$p #extracts p.values from dunn_test
names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
hold <-multcompLetters(holder) #this generates the letters 
letters.phiNPQ <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector
letters.phiNPQ

#create data frame for putting on letters
value_max_phiNPQ = d1 %>%
  group_by(Species) %>%
  summarize(max_value = max(phiNPQ))
value_max_phiNPQ

#create text frame for graphic overlay
text_phiNPQ <- cbind(value_max_phiNPQ, letters.phiNPQ) # combine both types of data 
text_phiNPQ=text_phiNPQ[,c(1,3,2)] #rearrange columns




#make plot
phiNPQ <- ggplot(d1, aes(x = Species, y = phiNPQ,
                         fill = Species)) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65, position = "jitter") +
  geom_violin(alpha = 0.3, width =1.25) +
  scale_y_continuous(limits = c(min(d1$phiNPQ)-.05, max(d1$phiNPQ)+.05)) +
  theme_bw() +
  theme(text = element_text(size = 24), axis.text.x = element_text(angle = 30, vjust = .6, hjust=.5),
        legend.position = "none", axis.title.x=element_blank(),
        legend.title = element_blank()) +
  scale_fill_manual(values = pal.2) +
  geom_text(data = text_phiNPQ, 
            aes(x=Species, y = 0.01 + max_value, label = value), 
            vjust=-0.3,
            size=7, 
            fontface = "bold") +
  labs(y=expression("Φ"[NPQ]))


#run nonparametric test: phiNO
(res.kruskal <- d1 %>% kruskal_test(phiNO ~ Species)) #this will indicate if there are siginficant differences between groups
d1 %>% kruskal_effsize(phiNO ~ Species) #this indicates the amount of variance between groups explained by the treatment. 
#The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
#value of .596 is very large it suggests that latitudinal sites vary by phiNO 

#Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
(pwc <- d1 %>% 
    dunn_test(phiNO ~ Species, p.adjust.method = "bonferroni")) #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 

#designate significant differences between groups with letters
holder <- pwc$p #extracts p.values from dunn_test
names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
hold <-multcompLetters(holder) #this generates the letters 
letters.phiNO <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector
letters.phiNO

#create data frame for putting on letters
value_max_phiNO = d1 %>%
  group_by(Species) %>%
  summarize(max_value = max(phiNO))
value_max_phiNO

#create text frame for graphic overlay
text_phiNO <- cbind(value_max_phiNO, letters.phiNO) # combine both types of data 
text_phiNO=text_phiNO[,c(1,3,2)] #rearrange columns




#make plot
phiNO <- ggplot(d1, aes(x = Species, y = phiNO,
                        fill = Species)) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65, position = "jitter") +
  geom_violin(alpha = 0.3, width =1.25) +
  scale_y_continuous(limits = c(min(d1$phiNO)-.05, max(d1$phiNO)+.05)) +
  theme_bw() +
  theme(text = element_text(size = 24), axis.text.x = element_text(angle = 30, vjust = .6, hjust=.5),
        legend.position = "none", axis.title.x=element_blank(),
        legend.title = element_blank()) +
  scale_fill_manual(values = pal.2) +
  geom_text(data = text_phiNO, 
            aes(x=Species, y = 0.01 + max_value, label = value), 
            vjust=-0.3,
            size=7, 
            fontface = "bold") +
  labs(y=expression("Φ"[NO]))

# #run nonparametric test: JETR
# (res.kruskal <- d1 %>% kruskal_test(JETR ~ Species)) #this will indicate if there are siginficant differences between groups
# d1 %>% kruskal_effsize(JETR ~ Species) #this indicates the amount of variance between groups explained by the treatment. 
# #The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
# #value of .704 is very large it suggests that latitudinal sites vary by JETR 
# 
# #Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
# (pwc <- d1 %>% 
#     dunn_test(JETR ~ Species, p.adjust.method = "bonferroni")) #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 
# 
# #designate significant differences between groups with letters
# holder <- pwc$p #extracts p.values from dunn_test
# names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
# hold <-multcompLetters(holder) #this generates the letters 
# letters.JETR <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector
# letters.JETR
# 
# #create data frame for putting on letters
# value_max_JETR = d1 %>%
#   group_by(Species) %>%
#   summarize(max_value = max(JETR))
# value_max_JETR
# 
# #create text frame for graphic overlay
# text_JETR <- cbind(value_max_JETR, letters.JETR) # combine both types of data 
# text_JETR=text_JETR[,c(1,3,2)] #rearrange columns
# 
# 
# 
# 
# #make plot
# JETR <- ggplot(d1, aes(x = Species, y = JETR,
#                        fill = Species)) +
#   geom_point(shape = 21, size = 3,
#              alpha = 0.65, position = "jitter") +
#   geom_violin(alpha = 0.3, width =1.25) +
#   scale_y_continuous(limits = c(min(d1$JETR)-.05, max(d1$JETR)+25)) +
#   theme_bw() +
#   theme(text = element_text(size = 24), axis.text.x = element_text(angle = 30, vjust = .6, hjust=.5),
#         legend.position = "none", axis.title.x=element_blank(),
#         legend.title = element_blank()) +
#   scale_fill_manual(values = pal.2) +
#   geom_text(data = text_JETR, 
#             aes(x=Species, y = 3 + max_value, label = value), 
#             vjust=-0.3,
#             size=7, 
#             fontface = "bold") +
#   labs(y=expression("J"[ETR]))
# 
# #run nonparametric test: JNPQ
# (res.kruskal <- d1 %>% kruskal_test(JNPQ ~ Species)) #this will indicate if there are siginficant differences between groups
# d1 %>% kruskal_effsize(JNPQ ~ Species) #this indicates the amount of variance between groups explained by the treatment. 
# #The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
# #value of .389 is very large it suggests that latitudinal sites vary by JNPQ 
# 
# #Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
# (pwc <- d1 %>% 
#     dunn_test(JNPQ ~ Species, p.adjust.method = "bonferroni")) #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 
# 
# #designate significant differences between groups with letters
# holder <- pwc$p #extracts p.values from dunn_test
# names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
# hold <-multcompLetters(holder) #this generates the letters 
# letters.JNPQ <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector
# letters.JNPQ
# 
# #create data frame for putting on letters
# value_max_JNPQ = d1 %>%
#   group_by(Species) %>%
#   summarize(max_value = max(JNPQ))
# value_max_JNPQ
# 
# #create text frame for graphic overlay
# text_JNPQ <- cbind(value_max_JNPQ, letters.JNPQ) # combine both types of data 
# text_JNPQ=text_JNPQ[,c(1,3,2)] #rearrange columns
# 
# 
# 
# 
# #make plot
# JNPQ <- ggplot(d1, aes(x = Species, y = JNPQ,
#                        fill = Species)) +
#   geom_point(shape = 21, size = 3,
#              alpha = 0.65, position = "jitter") +
#   geom_violin(alpha = 0.3, width =1.25) +
#   scale_y_continuous(limits = c(min(d1$JNPQ)-10, max(d1$JNPQ)+25)) +
#   theme_bw() +
#   theme(text = element_text(size = 24), axis.text.x = element_text(angle = 30, vjust = .6, hjust=.5),
#         legend.position = "none", axis.title.x=element_blank(),
#         legend.title = element_blank()) +
#   scale_fill_manual(values = pal.2) +
#   geom_text(data = text_JNPQ, 
#             aes(x=Species, y = 10 + max_value, label = value), 
#             vjust=-0.3,
#             size=7, 
#             fontface = "bold") +
#   labs(y=expression("J"[NPQ]))
# 
# #run nonparametric test: JNO
# (res.kruskal <- d1 %>% kruskal_test(JNO ~ Species)) #this will indicate if there are siginficant differences between groups
# d1 %>% kruskal_effsize(JNO ~ Species) #this indicates the amount of variance between groups explained by the treatment. 
# #The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
# #value of .677 is very large it suggests that latitudinal sites vary by JNO 
# 
# #Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
# (pwc <- d1 %>% 
#     dunn_test(JNO ~ Species, p.adjust.method = "bonferroni")) #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 
# 
# #designate significant differences between groups with letters
# holder <- pwc$p #extracts p.values from dunn_test
# names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
# hold <-multcompLetters(holder) #this generates the letters 
# letters.JNO <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector
# letters.JNO
# 
# #create data frame for putting on letters
# value_max_JNO = d1 %>%
#   group_by(Species) %>%
#   summarize(max_value = max(JNO))
# value_max_JNO
# 
# #create text frame for graphic overlay
# text_JNO <- cbind(value_max_JNO, letters.JNO) # combine both types of data 
# text_JNO=text_JNO[,c(1,3,2)] #rearrange columns
# 
# 
# 
# 
# #make plot
# JNO <- ggplot(d1, aes(x = Species, y = JNO,
#                       fill = Species)) +
#   geom_point(shape = 21, size = 3,
#              alpha = 0.65, position = "jitter") +
#   geom_violin(alpha = 0.3, width =1.25) +
#   scale_y_continuous(limits = c(min(d1$JNO)-10, max(d1$JNO)+25)) +
#   theme_bw() +
#   theme(text = element_text(size = 24), axis.text.x = element_text(angle = 30, vjust = .6, hjust=.5),
#         legend.position = "none", axis.title.x=element_blank(),
#         legend.title = element_blank()) +
#   scale_fill_manual(values = pal.2) +
#   geom_text(data = text_JNO, 
#             aes(x=Species, y = 10 + max_value, label = value), 
#             vjust=-0.3,
#             size=7, 
#             fontface = "bold") +
#   labs(y=expression("J"[NO]))
# 
# 
# d2 <- d1[!d1$Species=="Zea",]
# #fligner.test(FvFmprime ~ Species, center=median, d2)
# 
# #run nonparametric test: Jtot
# (res.kruskal <- d1 %>% kruskal_test(Jtot ~ Species)) #this will indicate if there are siginficant differences between groups
# d1 %>% kruskal_effsize(Jtot ~ Species) #this indicates the amount of variance between groups explained by the treatment. 
# #The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
# #value of .721 is very large it suggests that latitudinal sites vary by Jtot 
# 
# #Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
# (pwc <- d1 %>% 
#     dunn_test(Jtot ~ Species, p.adjust.method = "bonferroni")) #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 
# 
# #designate significant differences between groups with letters
# holder <- pwc$p #extracts p.values from dunn_test
# names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
# hold <-multcompLetters(holder) #this generates the letters 
# letters.Jtot <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector
# letters.Jtot
# 
# #create data frame for putting on letters
# value_max_Jtot = d1 %>%
#   group_by(Species) %>%
#   summarize(max_value = max(Jtot))
# value_max_Jtot
# 
# #create text frame for graphic overlay
# text_Jtot <- cbind(value_max_Jtot, letters.Jtot) # combine both types of data 
# text_Jtot=text_Jtot[,c(1,3,2)] #rearrange columns
# 
# 
# 
# 
# #make plot
# Jtot <- ggplot(d1, aes(x = Species, y = Jtot,
#                        fill = Species)) +
#   geom_point(shape = 21, size = 3,
#              alpha = 0.65, position = "jitter") +
#   geom_violin(alpha = 0.3, width =1.25) +
#   scale_y_continuous(limits = c(min(d1$Jtot)-1, max(d1$Jtot)+50)) +
#   theme_bw() +
#   theme(text = element_text(size = 24), axis.text.x = element_text(angle = 30, vjust = .6, hjust=.5),
#         legend.position = "none", axis.title.x=element_blank(),
#         legend.title = element_blank()) +
#   scale_fill_manual(values = pal.2) +
#   geom_text(data = text_Jtot, 
#             aes(x=Species, y = 15 + max_value, label = value), 
#             vjust=-0.3,
#             size=7, 
#             fontface = "bold") +
#   labs(y=expression("J"[tot]))

#make fake plot for legend
library(cowplot)
library(grid)
library(gridExtra) 
colnames(d1)[1] <- "Taxa"
colnames(text_FvFm_mean)[1] <- "Taxa"
filler <- ggplot(d1, aes(x = Taxa, y = FvFm_mean,
                         fill = Taxa)) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65, position = "jitter") +
  geom_violin(alpha = 0.3, width =1.25) +
  scale_y_continuous(limits = c(min(d1$FvFm_mean)-.05, max(d1$FvFm_mean)+.05)) +
  theme_bw() +
  scale_fill_manual(values = pal.2) +
  geom_text(data = text_FvFm_mean, 
            aes(x=Taxa, y = 0.02 + max_value, label = value), 
            vjust=-0.3,
            size=7, 
            fontface = "bold") +
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
EFGX'

# tiff("Fluorescence by taxa.tiff", units="in", width=20, height=10, res=300)
# wrap_plots(A = FvFm_mean, B = qP, C = qL, D= NPQ, E = phiPSII, F = phiNPQ,
#            G = phiNO, X = legend, design=layout) + plot_annotation(tag_levels = 'A')
# dev.off()

