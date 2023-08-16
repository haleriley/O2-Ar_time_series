## A more correct approach would be to derive a "true" O2 value for the MIMS
## based on existing calibrations, then base the model only on timepoints where
## the mims and miniDOT O2 values agree
getwd()
setwd('C://Users/haler/Documents/PhD-Bowman/O2-Ar_time_series/R_Data/')
set.seed(1234)

library.path <- .libPaths()[1]

# ---- library ----

library(tidyverse, lib.loc = library.path)
library(lubridate)
library(zoo)
library(aod)
library(ranger)
library(gbm)


col1.aou <- "#648fff"
col2.o2bio <- "#dc267f"
col3.aoucor <- "#fe6100" 
col4.aoupred <- "#785ef0" 
col5.other <- "#ffb000"


## Define a function for O2 saturation, mean salinity for the time-series is
## hard coded.

O2sat <- function(s, t){
  TS = log((298.15 - t) /  (273.15 + t))
  
  A0 = 5.80818
  A1 = 3.20684
  A2 = 4.11890
  A3 = 4.93845
  A4 = 1.01567
  A5 = 1.41575
  B0 = -7.01211 * 10 ** -3
  B1 = -7.25958 * 10 ** -3
  B2 = -7.93334 * 10 ** -3
  B3 = -5.54491 * 10 ** -3
  C0 = -1.32412 * 10 ** -7
  
  O2 = exp(A0 +
             A1*TS +
             A2* TS ** 2 +
             A3* TS ** 3 +
             A4* TS ** 4 +
             A5* TS ** 5 +
             s*(B0 +
                  B1* TS +
                  B2* TS ** 2 +
                  B3* TS ** 3 +
                  C0* s ** 2))
  
  return(O2)
  
}


# ---- train and validate boosted regression tree ----

combined.df <- readRDS(file = "2023-08-08_combined_env_data_hourly.rds")

predictors <- c('aop',
                #'delta',
                "temperature.sccoos",
                'salinity',
                "pressure",
                "WSPD",
                "GST",
                "WVHT",
                "DPD",
                #"ATMP",
                #"PRES",
                "o2_bio")

accessory <- c('Date.Time', 'Dissolved.Oxygen', 'O2')


## data selection for modeling, constrain by vacuum pressure (MIMS) ##
combined.df.select <- combined.df[which(combined.df$O2 >= 1.5e-9), c(predictors, accessory)]

## not sure why there are missing values, but seem to be the case!
combined.df.select <- na.omit(combined.df.select)



sqrt(mean((combined.df.select$o2_bio - combined.df.select$aop)^2))

random.dates <- sample(combined.df.select$Date.Time, size = 100, replace = F)

my.rmse.values <- as.data.frame(matrix(data = NA, nrow = 100, ncol = 6))
colnames(my.rmse.values) <- c("Date.Time", "RMSE_no_model", "RMSE_model", "RMSE_model_full", "RMSE_model_hypertuned", "RMSE_hyper_full")

my.rmse.values$Date.Time <- random.dates
my.rmse.values$RMSE_no_model <- sqrt(mean((combined.df.select$o2_bio - combined.df.select$aop)^2))

for(d in 1:length(random.dates)){
  
  ## date range for validation
  
  date1 <- random.dates[d]
  date2 <- random.dates[d] + 2*604800
  
  combined.df.train <- combined.df.select
  
  combined.df.test <- combined.df.select[which(combined.df.select$Date.Time > date1 &
                                                 combined.df.select$Date.Time < date2),]
  
  combined.df.train <- combined.df.train[which(!combined.df.train$Date.Time %in% combined.df.test$Date.Time),]
  
  combined.df.test.Date.Time <- combined.df.test$Date.Time
  combined.df.train.Date.Time <- combined.df.train$Date.Time
  
  combined.df.train <- combined.df.train[,predictors]
  combined.df.test <- combined.df.test[,predictors]
  
  
  
  # ---- train partial model ----
  
  combined.gbm <- gbm(o2_bio ~ .,
                      data = combined.df.train,
                      distribution = 'gaussian',
                      interaction.depth = 2,
                      n.trees = 4000,
                      bag.fraction = 0.5,
                      #train.fraction = 0.7,
                      shrinkage = 0.01,
                      verbose = T,
                      cv.folds = 2)
  
  print(combined.gbm)
  
  # ---- apply model to validation data ----
  
  aop.cor <- predict(combined.gbm, combined.df.test)
  
  aop.cor.delta <- aop.cor - combined.df.test$o2_bio #difference between predicted and o2bio

  aop.cor.rmse <- sqrt(mean(aop.cor.delta^2))
  
  my.rmse.values$RMSE_model[d] <- aop.cor.rmse
  
  full.time.series <- predict(combined.gbm, combined.df.select)
  full.rmse <- sqrt(mean((full.time.series - combined.df.select$o2_bio)^2))
  other.rmse <- sqrt(mean((full.time.series - combined.df.select$o2_bio)^2))
  
  my.rmse.values$RMSE_model_full[d] <- other.rmse
  
  
  # ---- parameter hypertuning ----
  
  hyper.grid <- expand.grid(
    interaction.depth = 2:6,
    n.trees = seq(1000, 4000, 1000),
    bag.fraction = seq(0.3, 0.7, 0.1),
    shrinkage = c(0.001, 0.01, 0.1),
    RMSE = NA
  )
  
  for(i in 1:nrow(hyper.grid)){
    
    ## try clause necessary because some parameter combinations are
    ## incompatible
    
    try({
      
      temp.model <- gbm(
        formula = o2_bio ~ .,
        data = combined.df.train,
        interaction.depth = hyper.grid$interaction.depth[i],
        distribution = 'gaussian',
        n.trees       = hyper.grid$n.trees[i],
        bag.fraction  = hyper.grid$bag.fraction[i],
        shrinkage   = hyper.grid$shrinkage[i],
        cv.folds = 2
      )
      
      temp.predict <- predict(temp.model, combined.df.test)
      temp.delta <- temp.predict - combined.df.test$o2_bio
      hyper.grid$RMSE[i] <- sqrt(mean(temp.delta^2))
    }, silent = F)
    
    print(paste(i, 'out of', nrow(hyper.grid), hyper.grid$RMSE[i]))
  }
  
  hyper.grid <- na.omit(hyper.grid)
  
  hist(hyper.grid$RMSE, breaks = 100)
  
  final.gbm <- gbm(
    formula = o2_bio ~ .,
    data = combined.df.train,
    interaction.depth = hyper.grid$interaction.depth[which.min(hyper.grid$RMSE)],
    distribution = 'gaussian',
    n.trees       = hyper.grid$n.trees[which.min(hyper.grid$RMSE)],
    bag.fraction  = hyper.grid$bag.fraction[which.min(hyper.grid$RMSE)],
    shrinkage   = hyper.grid$shrinkage[which.min(hyper.grid$RMSE)],
    cv.folds = 2
  )
  
  
  
  ## apply final model to test data
  
  aop.cor.final <- predict(final.gbm, combined.df.test)
  
  # ---- calculate model performance ----
  
  my.rmse <- sqrt(mean((combined.df.test$o2_bio - aop.cor.final)^2))
  # my.rmse.old <- sqrt(mean((combined.df.test$o2_bio - combined.df.test$aop)^2))
  
  my.rmse.values$RMSE_model_hypertuned[d] <- my.rmse
  
  
  # ---- full model ----
  
  full.gbm <- gbm(o2_bio ~ .,
                  data = na.omit(combined.df[,predictors]),
                  interaction.depth = hyper.grid$interaction.depth[which.min(hyper.grid$RMSE)],
                  distribution = 'gaussian',
                  n.trees       = hyper.grid$n.trees[which.min(hyper.grid$RMSE)],
                  bag.fraction  = hyper.grid$bag.fraction[which.min(hyper.grid$RMSE)],
                  shrinkage   = hyper.grid$shrinkage[which.min(hyper.grid$RMSE)],
                  cv.folds = 2)
  
  # save(list = c('combined.gbm', 'full.gbm'), file = '20230808_model.Rdata')
  
  
  
  # ---- apply model to miniDot timeseries ----
  
  full.predictors <- na.omit(combined.df.select[,c(full.gbm$var.names, 'Date.Time')])
  
  full.predictors$aop.corrected <- predict(full.gbm, full.predictors)
  
  my.rmse.values$RMSE_hyper_full[d] <- sqrt(mean((combined.df.select$o2_bio - full.predictors$aop.corrected)^2))
  
}

saveRDS(my.rmse.values, "2023-08-15_staggered_2week_RMSE_values.rds")

