#Kruskal-Wallace and Dunn testing for spectral data -- by taxa
#chris goodall -- 09.14.21
#re-wrote -- 01/24/22 
#updated species to taxa -- 10/09/2022
#updated taxa factor 04-19-2023

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
library(mosaic)
library(readr)

#import data
d1 <- read.csv("data/reflectance by taxa.csv")

#Clean up data to exclude physiologically impossible values
d1 <- d1[abs(d1$PRI) <= 1 &
           d1$ARI > -20 &
           d1$RGR > 0, ]

#remove sd columns
d1 <- d1[, -grep("_sd", colnames(d1))] 

d1$Taxa <- factor(d1$Taxa, levels = c("Sphagnum", "Marchantia", "Zea")) #re-order data so sphagnum is read first
#levels(d1$Taxa)

#run nonparametric test: NDVI
(res.kruskal <- d1 %>% kruskal_test(NDVI ~ Taxa)) #this will indicate if there are siginficant differences between groups
d1 %>% kruskal_effsize(NDVI ~ Taxa) #this indicates the amount of variance between groups explained by the treatment. 
#The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
#value of .553 is very large, it suggests that Taxa vary by NDVI 

#Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
(pwc <- d1 %>% 
    dunn_test(NDVI ~ Taxa, p.adjust.method = "bonferroni")) #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 
pwc

#designate significant differences between groups with letters
holder <- pwc$p #extracts p.values from dunn_test
names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
hold <-multcompLetters(holder) #this generates the letters 
letters.NDVI <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector

#create data frame for putting on letters
value_max_NDVI = d1 %>%
  group_by(Taxa) %>%
  summarize(max_value = max(NDVI))
value_max_NDVI

#create text frame for graphic overlay
text_NDVI <- cbind(value_max_NDVI, letters.NDVI) # combine both types of data 
text_NDVI=text_NDVI[,c(1,3,2)] #rearrange columns

#create color palette
pal <- plasma(3)[c(2,1,3)]

#make plot
NDVI <- ggplot(d1, aes(x = Taxa, y = NDVI,
                       fill = Taxa)) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65, position = "jitter") +
  geom_violin(alpha = 0.3, width =1.25) +
  scale_y_continuous(limits = c(min(d1$NDVI)-.1, max(d1$NDVI)+.1)) +
  theme_bw() +
  theme(text = element_text(size = 28), axis.text.x = element_text(angle = 30, vjust = .6, hjust=.5),
        legend.position = "none", axis.title.x=element_blank(),
        legend.title = element_blank()) +
  scale_fill_manual(values = pal) +
  geom_text(data = text_NDVI, 
            aes(x=Taxa, y = 0.02 + max_value, label = value), 
            vjust=-0.3,
            size=7, 
            fontface = "bold") 
NDVI


#run nonparametric test: PRI
(res.kruskal <- d1 %>% kruskal_test(PRI ~ Taxa)) #this will indicate if there are siginficant differences between groups
d1 %>% kruskal_effsize(PRI ~ Taxa) #this indicates the amount of variance between groups explained by the treatment. 
#The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
#value of .780 is very large, it suggests that Taxa vary by PRI 

#Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
(pwc <- d1 %>% 
    dunn_test(PRI ~ Taxa, p.adjust.method = "bonferroni")) #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 

#designate significant differences between groups with letters
holder <- pwc$p #extracts p.values from dunn_test
names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
hold <-multcompLetters(holder) #this generates the letters 
letters.PRI <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector

#create data frame for putting on letters
value_max_PRI = d1 %>%
  group_by(Taxa) %>%
  summarize(max_value = max(PRI))
value_max_PRI

#create text frame for graphic overlay
text_PRI <- cbind(value_max_PRI, letters.PRI) # combine both types of data 
text_PRI=text_PRI[,c(1,3,2)] #rearrange columns


#make plot
PRI <- ggplot(d1, aes(x = Taxa, y = PRI,
                      fill = Taxa)) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65, position = "jitter") +
  geom_violin(alpha = 0.3, width =1.25) +
  scale_y_continuous(limits = c(min(d1$PRI)-.1, max(d1$PRI)+.1)) +
  theme_bw() +
  theme(text = element_text(size= 28), axis.text.x = element_text(angle = 30, vjust = .6, hjust=.5),
        legend.position = "none", axis.title.x=element_blank(),
        legend.title = element_blank()) +
  scale_fill_manual(values = pal) +
  geom_text(data = text_PRI, 
            aes(x=Taxa, y = 0.02 + max_value, label = value), 
            vjust=-0.3,
            size=7, 
            fontface = "bold") 
PRI


#run nonparametric test: ExGm
(res.kruskal <- d1 %>% kruskal_test(ExGm ~ Taxa)) #this will indicate if there are siginficant differences between groups
d1 %>% kruskal_effsize(ExGm ~ Taxa) #this indicates the amount of variance between groups explained by the treatment. 
#The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
#value of .819 is very large, it suggests that Taxa vary by ExGm 

#Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
(pwc <- d1 %>% 
    dunn_test(ExGm ~ Taxa, p.adjust.method = "bonferroni")) #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 

#designate significant differences between groups with letters
holder <- pwc$p #extracts p.values from dunn_test
names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
hold <-multcompLetters(holder) #this generates the letters 
letters.ExGm <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector

#create data frame for putting on letters
value_max_PRI = d1 %>%
  group_by(Taxa) %>%
  summarize(max_value = max(ExGm))
value_max_PRI

#create text frame for graphic overlay
text_PRI <- cbind(value_max_PRI, letters.ExGm) # combine both types of data 
text_PRI=text_PRI[,c(1,3,2)] #rearrange columns

#make plot
ExGm <- ggplot(d1, aes(x = Taxa, y = ExGm,
                       fill = Taxa)) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65, position = "jitter") +
  geom_violin(alpha = 0.3, width =1.25) +
  scale_y_continuous(limits = c(min(d1$ExGm)-.1, max(d1$ExGm)+1)) +
  theme_bw() +
  theme(text = element_text(size= 28), axis.text.x = element_text(angle = 30, vjust = .6, hjust=.5),
        legend.position = "none", axis.title.x=element_blank(),
        legend.title = element_blank()) +
  scale_fill_manual(values = pal) +
  labs(y=expression('ExG'[M])) +
  geom_text(data = text_PRI, 
            aes(x=Taxa, y = 0.5 + max_value, label = value), 
            vjust=-0.3,
            size=7, 
            fontface = "bold") 
ExGm


#run nonparametric test: GRVI
(res.kruskal <- d1 %>% kruskal_test(GRVI ~ Taxa)) #this will indicate if there are siginficant differences between groups
d1 %>% kruskal_effsize(GRVI ~ Taxa) #this indicates the amount of variance between groups explained by the treatment. 
#The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
#value of .843 is very large, it suggests that Taxa vary by GRVI 

#Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
(pwc <- d1 %>% 
    dunn_test(GRVI ~ Taxa, p.adjust.method = "bonferroni")) #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 

#designate significant differences between groups with letters
holder <- pwc$p #extracts p.values from dunn_test
names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
hold <-multcompLetters(holder) #this generates the letters 
letters.GRVI <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector

#create data frame for putting on letters
value_max_GRVI = d1 %>%
  group_by(Taxa) %>%
  summarize(max_value = max(GRVI))
value_max_GRVI

#create text frame for graphic overlay
text_GRVI <- cbind(value_max_GRVI, letters.GRVI) # combine both types of data 
text_GRVI=text_GRVI[,c(1,3,2)] #rearrange columns

#make plot
GRVI <- ggplot(d1, aes(x = Taxa, y = GRVI,
                       fill = Taxa)) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65, position = "jitter") +
  geom_violin(alpha = 0.3, width =1.25) +
  scale_y_continuous(limits = c(min(d1$GRVI)-.1, max(d1$GRVI)+.5)) +
  theme_bw() +
  theme(text = element_text(size= 28), axis.text.x = element_text(angle = 30, vjust = .6, hjust=.5),
        legend.position = "none", axis.title.x=element_blank(),
        legend.title = element_blank()) +
  scale_fill_manual(values = pal) +
  geom_text(data = text_GRVI, 
            aes(x=Taxa, y = 0.15 + max_value, label = value), 
            vjust=-0.3,
            size=7, 
            fontface = "bold") 
GRVI


#run nonparametric test: ARI
(res.kruskal <- d1 %>% kruskal_test(ARI ~ Taxa)) #this will indicate if there are siginficant differences between groups
d1 %>% kruskal_effsize(ARI ~ Taxa) #this indicates the amount of variance between groups explained by the treatment. 
#The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
#value of .842 is very large, it suggests that Taxa vary by ARI 

#Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
(pwc <- d1 %>% 
    dunn_test(ARI ~ Taxa, p.adjust.method = "bonferroni")) #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 

#designate significant differences between groups with letters
holder <- pwc$p #extracts p.values from dunn_test
names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
hold <-multcompLetters(holder) #this generates the letters 
letters.ARI <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector

#create data frame for putting on letters
value_max_ARI = d1 %>%
  group_by(Taxa) %>%
  summarize(max_value = max(ARI))
value_max_ARI

#create text frame for graphic overlay
text_ARI <- cbind(value_max_ARI, letters.ARI) # combine both types of data 
text_ARI=text_ARI[,c(1,3,2)] #rearrange columns

#make plot
ARI <- ggplot(d1, aes(x = Taxa, y = ARI,
                      fill = Taxa)) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65, position = "jitter") +
  geom_violin(alpha = 0.3, width =1.25) +
  scale_y_continuous(limits = c(min(d1$ARI)-.1, max(d1$ARI)+.5)) +
  theme_bw() +
  theme(text = element_text(size= 28), axis.text.x = element_text(angle = 30, vjust = .6, hjust=.5),
        legend.position = "none", axis.title.x=element_blank(),
        legend.title = element_blank()) +
  scale_fill_manual(values = pal) +
  geom_text(data = text_ARI, 
            aes(x=Taxa, y = 0.15 + max_value, label = value), 
            vjust=-0.3,
            size=7, 
            fontface = "bold") 
ARI


#run nonparametric test: RGR
(res.kruskal <- d1 %>% kruskal_test(RGR ~ Taxa)) #this will indicate if there are siginficant differences between groups
d1 %>% kruskal_effsize(RGR ~ Taxa) #this indicates the amount of variance between groups explained by the treatment. 
#The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
#value of .846 is very large, it suggests that Taxa vary by RGR 

#Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
(pwc <- d1 %>% 
    dunn_test(RGR ~ Taxa, p.adjust.method = "bonferroni")) #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 

#designate significant differences between groups with letters
holder <- pwc$p #extracts p.values from dunn_test
names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
hold <-multcompLetters(holder) #this generates the letters 
letters.RGR <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector

#create data frame for putting on letters
value_max_RGR = d1 %>%
  group_by(Taxa) %>%
  summarize(max_value = max(RGR))
value_max_RGR

#create text frame for graphic overlay
text_RGR <- cbind(value_max_RGR, letters.RGR) # combine both types of data 
text_RGR=text_RGR[,c(1,3,2)] #rearrange columns


#make plot
RGR <- ggplot(d1, aes(x = Taxa, y = RGR,
                      fill = Taxa)) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65, position = "jitter") +
  geom_violin(alpha = 0.3, width =1.25) +
  scale_y_continuous(limits = c(min(d1$RGR)-.1, max(d1$RGR)+.5)) +
  theme_bw() +
  theme(text = element_text(size= 28), axis.text.x = element_text(angle = 30, vjust = .6, hjust=.5),
        legend.position = "none", axis.title.x=element_blank(),
        legend.title = element_blank()) +
  scale_fill_manual(values = pal) +
  geom_text(data = text_RGR, 
            aes(x=Taxa, y = 0.15 + max_value, label = value), 
            vjust=-0.3,
            size=7, 
            fontface = "bold") 
RGR

library(cowplot)
library(grid)
library(gridExtra) 
filler <- ggplot(d1, aes(x = Taxa, y = NDVI, #arbitrary index, doesn't matter
                         fill = Taxa)) +
  scale_fill_manual(values = pal) +
  #scale_y_continuous(limits = c(0, 50)) +
  geom_violin(alpha = 0.3, width=.6) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65,
             position = "jitter") +
  theme_bw() +
  theme(text = element_text(size= 28),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        legend.text=element_text(size=26),
        legend.title=element_text(size=32, face = "bold"),
        legend.position = "right",) 
legend <- cowplot::get_legend(filler)
filler
#make layout
library(patchwork)
layout <- '
abcg
abcg
abcg
defg
defg
defg
'
#  
# tiff("reflectance data by taxa.tiff", units="in", width=20, height=11.5, res=300)
# wrap_plots(a=NDVI,b=PRI,c=ARI, d=ExGm, e=GRVI, f=RGR, g=legend, design=layout) + plot_annotation(tag_levels = "A")
# dev.off()
