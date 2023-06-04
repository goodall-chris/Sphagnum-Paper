#Kruskal-Wallace and Dunn testing for spectral data -- latitude
#chris goodall -- 09.14.21
#re-wrote -- 01/24/22 
#updated 04-23-2023

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

#import data
d1 <- read.csv("reflectance by latitude 04-22-2023.csv")

#clean data
d1$lat <- round(d1$lat, 0) #round by lat
d1<- d1[, -which(names(d1) %in% c("X"))] #eliminate nonsense column

#Clean up data
d1 <- d1[abs(d1$PRI) <= 1 &
           d1$ARI > -20 &
           d1$RGR > 0, ]

colnames(d1)[2] <- "Latitude" #longer name, but makes graphics look better

#remove sd columns
d1 <- d1[, -grep("_sd", colnames(d1))] 


#run nonparametric test: NDVI
(res.kruskal <- d1 %>% kruskal_test(NDVI ~ Latitude)) #this will indicate if there are siginficant differences between groups
d1 %>% kruskal_effsize(NDVI ~ Latitude) #this indicates the amount of variance between groups explained by the treatment. 
#The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
#value of .579 is very large it suggests that latitudinal sites vary by NDVI 

#Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
pwc1<- (pwc <- d1 %>% 
    dunn_test(NDVI ~ Latitude, p.adjust.method = "bonferroni")) #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 

#designate significant differences between groups with letters
holder <- pwc$p #extracts p.values from dunn_test
names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
hold <-multcompLetters(holder) #this generates the letters 
letters.NDVI <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector

#create data frame for putting on letters
value_max_NDVI = d1 %>%
  group_by(Latitude) %>%
  summarize(max_value = max(NDVI))
value_max_NDVI

#create text frame for graphic overlay
text_NDVI <- cbind(value_max_NDVI, letters.NDVI) # combine both types of data 
text_NDVI=text_NDVI[,c(1,3,2)] #rearrange columns

#add sneaky features for plotting
d1$lat <- as.factor(d1$Latitude)
text_NDVI$lat <- as.factor(text_NDVI$Latitude) 

#color palette
pal.2 <- viridis(9, direction = 1) #create palette

#make plot
NDVI <- ggplot(d1, aes(x = Latitude, y = NDVI,
                       fill = lat)) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65, position = "jitter") +
  geom_violin(alpha = 0.3, width =1.25) +
  scale_y_continuous(limits = c(min(d1$NDVI)-.1, max(d1$NDVI)+.1)) +
  theme_bw() +
  theme(text = element_text(size = 24),
        legend.position = "none",
        legend.title = element_blank()) +
  scale_fill_manual(values = pal.2) +
  geom_text(data = text_NDVI, 
            aes(x=Latitude, y = 0.02 + max_value, label = value), 
            vjust=-0.3,
            size=5.5, 
            fontface = "bold") 
NDVI

#run nonparametric test: PRI
(res.kruskal <- d1 %>% kruskal_test(PRI ~ Latitude)) #this will indicate if there are siginficant differences between groups
d1 %>% kruskal_effsize(PRI ~ Latitude) #this indicates the amount of variance between groups explained by the treatment. 
#The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
#value of .194 is large it suggests that latitudinal sites vary by PRI 

#Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
(pwc <- d1 %>% 
    dunn_test(PRI ~ Latitude, p.adjust.method = "bonferroni")) #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 

#designate significant differences between groups with letters
holder <- pwc$p #extracts p.values from dunn_test
names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
hold <-multcompLetters(holder) #this generates the letters 
letters.PRI <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector

#create data frame for putting on letters
value_max_PRI = d1 %>%
  group_by(Latitude) %>%
  summarize(max_value = max(PRI))
value_max_PRI

#create text frame for graphic overlay
text_PRI <- cbind(value_max_PRI, letters.PRI) # combine both types of data 
text_PRI=text_PRI[,c(1,3,2)] #rearrange columns

#add sneaky features for plotting
d1$lat <- as.factor(d1$Latitude)
text_PRI$lat <- as.factor(text_PRI$Latitude) 

#adjust text to fix cluttering 
text_PRI.1 <- text_PRI[!text_PRI$Latitude==31,]
text_PRI.2 <- text_PRI[text_PRI$Latitude==31,]

#make plot
PRI <- ggplot(d1, aes(x = Latitude, y = PRI,
                      fill = lat)) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65, position = "jitter") +
  geom_violin(alpha = 0.3, width =1.25) +
  scale_y_continuous(limits = c(min(d1$PRI)-.1, max(d1$PRI)+.1)) +
  theme_bw() +
  theme(text = element_text(size = 24),
        legend.position = "none",
        legend.title = element_blank()) +
  scale_fill_manual(values = pal.2) +
  geom_text(data = text_PRI.1, 
            aes(x=Latitude, y = 0.02 + max_value, label = value), 
            vjust=-0.3,
            size=5.5, 
            fontface = "bold") +
 geom_text(data = text_PRI.2, 
            aes(x=Latitude, y = 0.04 + max_value, label = value), 
            vjust=-0.3,
            size=5.5, 
            fontface = "bold")
PRI

#run nonparametric test: ExGm
# --by individual--
(res.kruskal <- d1 %>% kruskal_test(ExGm ~ Latitude)) #this will indicate if there are siginficant differences between groups
d1 %>% kruskal_effsize(ExGm ~ Latitude) #this indicates the amount of variance between groups explained by the treatment. 
#The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
#value of .396 is very large it suggests that latitudinal sites vary by ExGm 

#Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
(pwc <- d1 %>% 
    dunn_test(ExGm ~ Latitude, p.adjust.method = "bonferroni")) #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 

#designate significant differences between groups with letters
holder <- pwc$p #extracts p.values from dunn_test
names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
hold <-multcompLetters(holder) #this generates the letters 
letters.ExGm <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector

#create data frame for putting on letters
value_max_ExGm = d1 %>%
  group_by(Latitude) %>%
  summarize(max_value = max(ExGm))
value_max_ExGm

#create text frame for graphic overlay
text_ExGm <- cbind(value_max_ExGm, letters.ExGm) # combine both types of data 
text_ExGm=text_ExGm[,c(1,3,2)] #rearrange columns

#add sneaky features for plotting
d1$lat <- as.factor(d1$Latitude)
text_ExGm$lat <- as.factor(text_ExGm$Latitude) 

#adjust text to fix cluttering 
text_ExGm.1 <- text_ExGm[!text_ExGm$Latitude==32,]
text_ExGm.2 <- text_ExGm[text_ExGm$Latitude==32,]

#make plot
ExGm <- ggplot(d1, aes(x = Latitude, y = ExGm,
                       fill = lat)) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65, position = "jitter") +
  geom_violin(alpha = 0.3, width =1.25) +
  scale_y_continuous(limits = c(min(d1$ExGm)-.6, max(d1$ExGm)+1)) +
  theme_bw() +
  theme(text = element_text(size = 24),
        legend.position = "none",
        legend.title = element_blank()) +
  scale_fill_manual(values = pal.2) +
  geom_text(data = text_ExGm.1, 
            aes(x=Latitude, y = 0.3 + max_value, label = value, ), 
            vjust=-.5,
            size=5.5, 
            fontface = "bold") +
  geom_text(data = text_ExGm.2, 
            aes(x=Latitude, y = 0.6 + max_value, label = value), 
            vjust=-.3,
            size=5.5, 
            fontface = "bold") 
ExGm


#run nonparametric test: GRVI
# --by individual--
(res.kruskal <- d1 %>% kruskal_test(GRVI ~ Latitude)) #this will indicate if there are siginficant differences between groups
d1 %>% kruskal_effsize(GRVI ~ Latitude) #this indicates the amount of variance between groups explained by the treatment. 
#The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
#value of .389 is very large it suggests that latitudinal sites vary by GRVI 

#Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
(pwc <- d1 %>% 
    dunn_test(GRVI ~ Latitude, p.adjust.method = "bonferroni")) #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 

#designate significant differences between groups with letters
holder <- pwc$p #extracts p.values from dunn_test
names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
hold <-multcompLetters(holder) #this generates the letters 
letters.GRVI <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector

#create data frame for putting on letters
value_max_GRVI = d1 %>%
  group_by(Latitude) %>%
  summarize(max_value = max(GRVI))
value_max_GRVI

#create text frame for graphic overlay
text_GRVI <- cbind(value_max_GRVI, letters.GRVI) # combine both types of data 
text_GRVI=text_GRVI[,c(1,3,2)] #rearrange columns

#add sneaky features for plotting
d1$lat <- as.factor(d1$Latitude)
text_GRVI$lat <- as.factor(text_GRVI$Latitude) 

#adjust text to fix cluttering 
text_GRVI.1 <- text_GRVI[!text_GRVI$Latitude==31,]
text_GRVI.2 <- text_GRVI[text_GRVI$Latitude==31,]

#make plot
GRVI <- ggplot(d1, aes(x = Latitude, y = GRVI,
                       fill = lat)) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65, position = "jitter") +
  geom_violin(alpha = 0.3, width = 1.25) +
  scale_y_continuous(limits = c(min(d1$GRVI)-.25, max(d1$GRVI)+.25)) +
  theme_bw() +
  theme(text = element_text(size = 24),
        legend.position = "none",
        legend.title = element_blank()) +
  scale_fill_manual(values = pal.2) +
  geom_text(data = text_GRVI.1, 
            aes(x=Latitude, y = 0.05 + max_value, label = value), 
            vjust=-0.3,
            size=5.5, 
            fontface = "bold")+
  geom_text(data = text_GRVI.2, 
            aes(x=Latitude, y = 0.1 + max_value, label = value), 
            vjust=-0.3,
            size=5.5, 
            fontface = "bold") 
GRVI


#run nonparametric test: ARI
# --by individual--
(res.kruskal <- d1 %>% kruskal_test(ARI ~ Latitude)) #this will indicate if there are siginficant differences between groups
d1 %>% kruskal_effsize(ARI ~ Latitude) #this indicates the amount of variance between groups explained by the treatment. 
#The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
#value of .454 is very large it suggests that latitudinal sites vary by ARI 

#Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
(pwc <- d1 %>% 
    dunn_test(ARI ~ Latitude, p.adjust.method = "bonferroni")) #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 

#designate significant differences between groups with letters
holder <- pwc$p #extracts p.values from dunn_test
names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
hold <-multcompLetters(holder) #this generates the letters 
letters.ARI <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector

#create data frame for putting on letters
value_max_ARI = d1 %>%
  group_by(Latitude) %>%
  summarize(max_value = max(ARI))
value_max_ARI

#create text frame for graphic overlay
text_ARI <- cbind(value_max_ARI, letters.ARI) # combine both types of data 
text_ARI=text_ARI[,c(1,3,2)] #rearrange columns

#add sneaky features for plotting
d1$lat <- as.factor(d1$Latitude)
text_ARI$lat <- as.factor(text_ARI$Latitude) 

#adjust text to fix cluttering 
text_ARI.1 <- text_ARI[!text_ARI$Latitude==31,]
text_ARI.2 <- text_ARI[text_ARI$Latitude==31,]

#make plot
ARI <- ggplot(d1, aes(x = Latitude, y = ARI,
                      fill = lat)) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65, position = "jitter") +
  geom_violin(alpha = 0.3, width = 1.25) +
  scale_y_continuous(limits = c(min(d1$ARI)-.25, max(d1$ARI)+.25)) +
  theme_bw() +
  theme(text = element_text(size = 24),
        legend.position = "none",
        legend.title = element_blank()) +
  scale_fill_manual(values = pal.2) +
  geom_text(data = text_ARI.1, 
            aes(x=Latitude, y = 0.05 + max_value, label = value), 
            vjust=-0.3,
            size=5.5, 
            fontface = "bold") +
  geom_text(data = text_ARI.2, 
            aes(x=Latitude, y = 0.05 + max_value, label = value), 
            vjust=-0.3,
            hjust=.75,
            size=5.5, 
            fontface = "bold") 
ARI


#run nonparametric test: RGR
# --by individual--
(res.kruskal <- d1 %>% kruskal_test(RGR ~ Latitude)) #this will indicate if there are siginficant differences between groups
d1 %>% kruskal_effsize(RGR ~ Latitude) #this indicates the amount of variance between groups explained by the treatment. 
#The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
#value of .412 is very large it suggests that latitudinal sites vary by RGR 

#Post hoc analysis: In this context, we want to use# a Dunn's test as it uses the same underlying assumptions of the KW test
(pwc <- d1 %>% 
    dunn_test(RGR ~ Latitude, p.adjust.method = "bonferroni")) #you can also set p.adjust.method = "fdr" here, but "bonferroni" involves more tight confidence intervals. 

#designate significant differences between groups with letters
holder <- pwc$p #extracts p.values from dunn_test
names(holder)  <-  paste(pwc$group1, pwc$group2, sep="-")
hold <-multcompLetters(holder) #this generates the letters 
letters.RGR <- reshape::melt(hold$Letters) #this pulls the letters from the list into a vector

#create data frame for putting on letters
value_max_RGR = d1 %>%
  group_by(Latitude) %>%
  summarize(max_value = max(RGR))
value_max_RGR

#create text frame for graphic overlay
text_RGR <- cbind(value_max_RGR, letters.RGR) # combine both types of data 
text_RGR=text_RGR[,c(1,3,2)] #rearrange columns

#add sneaky features for plotting
d1$lat <- as.factor(d1$Latitude)
text_RGR$lat <- as.factor(text_RGR$Latitude) 


#deal with overlap 
text_RGR.1 <- text_RGR[!text_RGR$Latitude==32,]
text_RGR.2 <- text_RGR[text_RGR$Latitude==32,]

#make plot
RGR <- ggplot(d1, aes(x = Latitude, y = RGR,
                      fill = lat)) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65, position = "jitter") +
  geom_violin(alpha = 0.3, width =2) +
  scale_y_continuous(limits = c(min(d1$RGR)-.25, max(d1$RGR)+.8)) +
  theme_bw() +
  theme(text = element_text(size = 24),
        legend.position = "none",
        legend.title = element_blank()) +
  scale_fill_manual(values = pal.2) +
  geom_text(data = text_RGR.1, 
                  aes(x=Latitude, y = max_value +.1, label = value), 
                  vjust=-0.3,
                  size=5.5, 
                  fontface = "bold") +
  geom_text(data = text_RGR.2, 
            aes(x=Latitude, y = max_value +.4, label = value), 
            vjust=-0.3,
            size=5.5, 
            fontface = "bold")


#make separatae figure for legend
d2 <- d1 
colnames(d2) <- c("population", "lat", "NDVI", "PRI", "ExGm", "GRVI", "ARI", "RGR", "Latitude")

library(cowplot)
library(grid)
library(gridExtra) 
filler <- ggplot(d2, aes(x = lat, y = NDVI, #arbitrary index, doesn't matter
                         fill = Latitude)) +
  scale_fill_manual(values = pal.2) +
  #scale_y_continuous(limits = c(0, 50)) +
  geom_violin(alpha = 0.3, width=.6) +
  geom_point(shape = 21, size = 3,
             alpha = 0.65,
             position = "jitter") +
  theme_bw() +
  theme(text = element_text(size = 24),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        legend.text=element_text(size=26),
        legend.title=element_text(size=32, face="bold"),
        legend.position = "right",) 
legend <- cowplot::get_legend(filler)

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

tiff("reflectance data by latitude -- 04-22-2023.tiff", units="in", width=20, height=10, res=300)
wrap_plots(a=NDVI,b=PRI,c=ExGm, d=GRVI, e=ARI, f=RGR, g=legend, design=layout) + plot_annotation(tag_levels = "A")
dev.off()

#insert ggplot code

