#FDA for reflectance -- Taxa
#Flexible Discriminant Analysis update -- 9.29.21
#cleaned 01/27/21
#updated taxa factor order 04-22-2023

#packages 
library(progress)
library(dplyr)
library(caret)
library(mda)
library(klaR)
library(sda)
library(viridis)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(readr)

#data 
d1 <- read_csv("data/reflectance by taxa 04-22-2023.csv")

#Clean up data to exclude physiologically impossible values
d1 <- d1[abs(d1$PRI) <= 1 &
           d1$ARI > -20 &
           d1$RGR > 0, ]

#remove sd columns
d1 <- d1[, -grep("_sd", colnames(d1))] 

d1$Taxa <- factor(d1$Taxa, levels = c("Sphagnum", "Marchantia", "Zea")) #re-order data so sphagnum is read first
#levels(d1$Taxa)

#creating data training partitions
trains <- c(0.4, 0.45, 0.5, 0.55, 0.6,
            0.65, 0.7, 0.75, 0.8)

#build ind level FDA by lat
pred_success_1 <- vector("list", 0)
pred_cond_1 <- vector("list", 0)
ival_1 <- vector("list", 0)

#progress bar
pb <- progress_bar$new(
  format = " running [:bar] :percent eta: :eta",
  total = 100, clear = FALSE, width= 60) #set up progress bar to put nerves to rest
#function
for(i in 1:100){
  set.seed(i)
  pb$tick()
  for(j in 1:length(trains)){
    data_train <- d1$Taxa %>%
      createDataPartition(p = trains[j], list = FALSE)#65%
    train.data <- d1 [data_train, ]
    test.data <- d1 [-data_train, ]
    preproc.param <- train.data %>%
      preProcess(method = c("center", "scale"))
    # Transform the data using the estimated parameters
    train.transformed <- preproc.param %>% predict(train.data)
    test.transformed <- preproc.param %>% predict(test.data)
    # Fit the model
    model <- fda(Taxa~. , data = train.transformed[, c(1:7)])
    # Make predictions
    predictions <- model %>% predict(test.transformed)
    # Model accuracy
    pred_success_1[[(j + (i - 1) * length(trains))]] <- mean(predictions==test.transformed$Taxa)
    pred_cond_1[[(j + (i - 1) * length(trains))]] <- trains[j]
    ival_1[[(j + (i - 1) * length(trains))]] <- i
    #mean(predictions==test.transformed$Taxa)
  }
}

pred_success_1 <- do.call("rbind", pred_success_1)
pred_cond_1 <- do.call("rbind", pred_cond_1)
ival_1 <- do.call("rbind", ival_1)
success <- data.frame(pred_success_1, pred_cond_1, ival_1)
summary(success)

#create data partitions for modeling using individual data and by Taxa
set.seed(26)
data_train <- d1$Taxa %>%
  createDataPartition(p = 0.6, list = FALSE) #60/40 split
train.data <- d1 [data_train, ]
test.data <- d1 [-data_train, ]
preproc.param <- train.data %>% 
  preProcess(method = c("center", "scale")) #calculates predictors

# Transform the data using the estimated parameters
train.transformed <- preproc.param %>% predict(train.data) #
test.transformed <- preproc.param %>% predict(test.data)

# Fit the model
model_1 <- fda(Taxa~. , data = train.transformed[, c(1:7)])
#see linear discriminant loadings for % variance explained
model_1
# Make predictions
predictions_1 <- model_1 %>% predict(test.transformed)

# Model accuracy
(accuracy_1<- mean(predictions_1==test.transformed$Taxa))
accuracy_1<- round(accuracy_1, 3)

#make data to plot 
plot.data_1 <- model_1$fit$fitted.values %>% 
  as_tibble() %>% 
  bind_cols(Taxa = train.transformed[,"Taxa"]) 

#function for formatting axis labels
percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}

#create plot 
pal <- plasma(3)[c(2,1,3)] #number of Taxa spots
figa <- ggplot(plot.data_1, aes(V1, V2)) +
  labs(x = paste("FDA1 (", percent(model_1$percent.explained[1]/100), ") \n Model Accuracy = ", percent(accuracy_1), sep=""),
       y = paste("FDA2 (", percent(model_1$percent.explained[2]/100 - model_1$percent.explained[1]/100), ")", sep=""),
       #title = "  Individual",
       fill = "Taxa",
       color = "Taxa") +
  ylim(-1.5,3)+ #change appropriately
  stat_density_2d(aes(fill = Taxa), alpha = 0.3, geom = "polygon", show.legend = T) +
  stat_ellipse(aes(x = V1, y = V2, color = Taxa),
               size = 1.5, show.legend = FALSE,
               type = "t",
               level = 0.9) +
  geom_point(aes(fill = Taxa, shape = popID),
             size = 3,
             shape = 21,
             alpha = 0.5) +
  #guides(fill = guide_legend(nrow = 5, byrow = TRUE, title="Taxa")) +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +
  theme_classic() +
  theme(text = element_text(size = 24),
        legend.position = "right",
        legend.title = element_text(face="bold")
  )
figa

tiff("FDA for reflectance data -- Taxa -- 04-22-2023.tiff", units="in", width=10, height=9, res=300)
figa
dev.off()
