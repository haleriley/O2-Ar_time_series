# Analysis of Environmental Conditions
# RJH
# 2022-02-03 

# ---- library ----

library(tidyverse)
library(lubridate)
library(zoo)

# ---- read in data ----

setwd("C://Users/haler/Documents/PhD-Bowman/MIMS-miniDOT_O2-Ar_Study/MIMS_4riley/")

mims.zoo <- readRDS(file = "2023-02-03_mims_zoo.rds")
miniDot.zoo <- readRDS(file = "2023-02-03_miniDot_zoo.rds")
sccoos.zoo <- readRDS(file = "2023-02-03_sccoos_zoo.rds")
ljac.zoo <- readRDS(file = "2023-02-03_ljac_zoo.rds")
ljpc1.zoo <- readRDS(file = "2023-02-03_ljpc1_zoo.rds")
# combined.df <- readRDS(file = "2022-10-10_combined_env_data_hourly.rds")
combined.df <- readRDS(file = "2023-02-03_combined.df.select.rds") # this is the dataset used in the model
aou.df <- readRDS(file = "full_aou_df.rds")

# ---- format zoo to data frames ----

mims.df <- data.frame(mims.zoo)
miniDot.df <- data.frame(miniDot.zoo)
sccoos.df <- data.frame(sccoos.zoo)
ljac.df <- data.frame(ljac.zoo)
ljpc1.df <- data.frame(ljpc1.zoo)

# ---- compare hourly data used in model to total hourly data ----

# how well does the model data represent the entire time series?
compare.model.to.all.data <- function(model.data, all.data){
  
  my.model.data <- as.numeric(model.data[which(is.na(model.data) == FALSE)])
  my.all.data <- as.numeric(all.data[which(is.na(all.data) == FALSE)])
  
  model.mean <- mean(my.model.data)
  model.sd <- sd(my.model.data)
  model.norm <- shapiro.test(my.model.data)
  all.mean <- mean(my.all.data)
  all.sd <- sd(my.all.data)
  all.norm <- shapiro.test(sample(my.all.data, size = 5000))
  
  my.t <- t.test(x = my.model.data, y = my.all.data, alternative = "two.sided")
  
  my.summary <- cbind(c("model.mean", "model.sd", "model.norm.p.value", "all.mean", "all.sd", "all.norm.p.value", "t.test.p.value"), c(model.mean, model.sd, model.norm$p.value, all.mean, all.sd, all.norm$p.value, my.t$p.value))
  
  curve(expr = dnorm(x, mean = model.mean, sd = model.sd), col = "blue", from = min((model.mean-(4*model.sd)), (all.mean-(4*all.sd))), to = max((model.mean+(4*model.sd)), (all.mean+(4*all.sd))), ylab = "Frequency")
  curve(expr = dnorm(x, mean = all.mean, sd = all.sd), col = "red", add = TRUE)
  legend("topright", legend = c("model", "all"), col = c("blue", "red"), lwd = 2)
  
  return(my.summary)
  
}

# for every two-week validation date range within the model data, how often does the validation data represent the model data (p > 0.05 t test)
compare.two.week.validation.data.to.model.data <- function(model.df = combined.df, model.col.name){
  
  model.df$dates <- as.Date(model.df$dates)
  unique.dates <- unique(combined.df$dates)
  
  my.t.values <- c()
  for(d in 1:(length(unique.dates))){
    
    my.date1 <- unique.dates[d]-6
    my.date2 <- unique.dates[d]+7
    
    my.validation.df <- model.df[which(model.df$dates >= my.date1 & model.df$dates <= my.date2),]
    my.validation.data <- my.validation.df[,which(colnames(my.validation.df) == model.col.name)]
    
    if(length(my.validation.data) > 1){
      
      my.model.data <- model.df[,which(colnames(model.df) == model.col.name)]
      
      my.t <- t.test(x = my.validation.data, y = my.model.data)
      
      my.t.values <- c(my.t.values, my.t$p.value)
      
    }
    
  }
  
  my.ratio <- length(which((my.t.values >= 0.05) == TRUE))/length(my.t.values)
  
  return(my.ratio)
  
}

find.trends.in.environmental.conditions <- function(full.data.df, full.data.colname){
  
  my.date <- parse_date_time(rownames(full.data.df), orders = "Ymd HMS")
  my.data <- as.data.frame(full.data.df[,which(colnames(full.data.df) == full.data.colname)])
  
  my.data <- data.frame(cbind(my.date, my.data))
  
  colnames(my.data) <- c("Date", full.data.colname)
  
  my.lm <- lm(my.data[,2]~my.data[,1])
  my.summary <- summary(my.lm)
  
  plot(my.data[,2]~my.data[,1], xlab = "Date", ylab = full.data.colname)
  abline(a = my.lm$coefficients[1], b = my.lm$coefficients[2], col = "red")
  
  return(my.summary)
  
}


# sccoos water temp
compare.model.to.all.data(model.data = combined.df$temperature.sccoos.zoo, all.data = sccoos.df$temperature)
# compare.two.week.validation.data.to.model.data(model.col.name = "temperature.sccoos.zoo")
find.trends.in.environmental.conditions(full.data.df = sccoos.df, full.data.colname = "temperature")

t.test(x = as.numeric(combined.df$temperature.sccoos.zoo), mu = mean(as.numeric(sccoos.df$temperature)), alternative = "two.sided")


# SCCOOS salinity
compare.model.to.all.data(combined.df$salinity, sccoos.df$salinity[which(sccoos.df$salinity > 33)]) # plot looks funky but that might be because we cherry picked salinity
# compare.two.week.validation.data.to.model.data(model.col.name = "salinity")
find.trends.in.environmental.conditions(full.data.df = sccoos.df[which(sccoos.df$salinity > 33),], full.data.colname = "salinity")


# SCCOOS pressure
compare.model.to.all.data(combined.df$pressure, sccoos.df$pressure)
# compare.two.week.validation.data.to.model.data(model.col.name = "pressure")
find.trends.in.environmental.conditions(full.data.df = sccoos.df, full.data.colname = "pressure")

# AOU
compare.model.to.all.data(combined.df$aou, aou.df$aou)
# find.trends.in.environmental.conditions(full.data.df = aou.df, full.data.colname = "aou")



# LJAC wind speed
compare.model.to.all.data(combined.df$WSPD, ljac.df$WSPD)
# compare.two.week.validation.data.to.model.data(model.col.name = "WSPD")
find.trends.in.environmental.conditions(full.data.df = ljac.df, full.data.colname = "WSPD")

# LJAC gust speed
compare.model.to.all.data(combined.df$GST, ljac.df$GST)
find.trends.in.environmental.conditions(full.data.df = ljac.df, full.data.colname = "GST")


my.colnames <- colnames(combined.df)[-10]

my.results <- cbind(my.colnames, rep(NA, times = length(my.colnames)))
for(c in 1:length(my.colnames)){
  
  my.name <- my.colnames[c]
  
  save.this <- compare.two.week.validation.data.to.model.data(model.col.name = my.name)
  
  my.results[c,2] <- save.this*100
  
}
val.data.summary <- my.results
# if you take a random two-week validation period within the model data, what is the probability that the validation data will have the same mean as the rest of the model data? 


