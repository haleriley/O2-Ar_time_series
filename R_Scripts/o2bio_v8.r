## A more correct approach would be to derive a "true" O2 value for the MIMS
## based on existing calibrations, then base the model only on timepoints where
## the mims and miniDOT O2 values agree
getwd()
setwd('C://Users/haler/Documents/PhD-Bowman/MIMS-miniDOT_O2-Ar_Study/MIMS_4riley/')
set.seed(1234)

## Define a function for O2 saturation, mean salinity for the time-series is
## hard coded.

O2sat <- function(s = 33.24, t){
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

#### get data ####

## Get O2bio calculated in read_lvm.py

mims <- read.csv('o2bio.csv') # this file is produced by read_lvm.py, refer to that script for details
mims$date_time <- strptime(mims$date_time, format = '%Y-%m-%d %H')
mims1 <- mims[mims$N2.Ar < 40 & mims$N2.Ar > 30 & mims$date_time < strptime('2021-3-26', format = '%Y-%m-%d'),]
mims2 <- mims[mims$N2.Ar < 20 & mims$N2.Ar > 9 & mims$date_time >= strptime('2021-3-26', format = '%Y-%m-%d'),]
mims <- rbind(mims1, mims2)
mims.date <- mims$date_time
mims$date_time <- NULL

mims.hourly <- aggregate(mims, by = list(as.character(mims.date)), FUN = mean, na.rm = T)
mims.hourly.date <- strptime(mims.hourly$Group.1, format = '%Y-%m-%d %H')

plot(mims.date, mims$o2_bio)
points(mims.hourly.date, mims.hourly$o2_bio,
       type = 'l',
       col = 'red')

plot(mims.date, mims$N2.Ar)

### Get miniDOT data ###

miniDot <- read.csv("SIOPier_miniDOT_20180809_20220920_JBowman.csv", header = F, na.strings = "NA", row.names = 1, skip = 2)
colnames(miniDot) <- c("UTC_Date_._Time", "Pacific.Standard.Time", "Temperature", "Dissolved.Oxygen", "Dissolved.Oxygen.Saturation", "Sensor")
# miniDot.colo <- read.csv('../SIOPier_miniDOT_20210430_20210709_BowmanLabSensor.csv', header = T, na.string = "NA")

miniDot$Pacific.Standard.Time <- strptime(miniDot$Pacific.Standard.Time, format = '%Y-%m-%d %H')

## MiniDot reports in mg/L, need uMol

miniDot$Dissolved.Oxygen <- (miniDot$Dissolved.Oxygen / (15.999 * 2 * 1000)) * 10 ** 6

miniDot.col.select <- c("Dissolved.Oxygen", "Temperature")

miniDot.hourly <- aggregate(miniDot[,miniDot.col.select], by = list(as.character(miniDot$Pacific.Standard.Time)), FUN = mean, na.rm = T)

### sccoos data ###

## this command should work but doesn't, download with wget instead
# sccoos <- read.csv('https://erddap.sccoos.org/erddap/tabledap/autoss.csv?station%2Ctime%2Ctemperature%2Cconductivity%2Cpressure%2Cchlorophyll%2Csalinity&station=%22scripps_pier%22')
# sccoos <- read.csv(file = "autoss_cd8b_3ba9_f9f5.csv")
# sccoos <- sccoos[-1,]
sccoos.columns <- c('station', 'time', 'temperature', 'conductivity','pressure','chlorophyll','salinity')
sccoos <- read.csv('sccoos.csv', header = T, skip = 2, col.names = sccoos.columns)
sccoos$station <- NULL
sccoos$time <- strptime(sccoos$time, format = '%Y-%m-%dT%H', tz = 'GMT')
sccoos$time <- as.POSIXlt(sccoos$time, tz = 'PST')
sccoos.col.select <- c('temperature', 'pressure','chlorophyll','salinity')
sccoos.hourly <- aggregate(sccoos[,sccoos.col.select], by = list(as.character(sccoos$time)), FUN = mean, na.rm = T)

### NOAA - ljac ###

## download historic data files with wget from https://www.ndbc.noaa.gov/station_history.php?station=ljac1

noaa.columns <- c('YY','MM','DD','hh','mm','WDIR','WSPD','GST','WVHT','DPD','APD','MWD','PRES','ATMP','WTMP','DEWP','VIS','TIDE')
noaa.na <- c('99.00', '999', '999.0', '9999.0', '99.0')

ljac.col.select <- c('WSPD','GST','PRES','ATMP','WTMP')
ljac <- read.table('NOAA_ljac/combined_ljac.txt', col.names = noaa.columns, na.strings = noaa.na)
ljac$date.time <- paste0(ljac$YY, '-', ljac$MM, '-', ljac$DD, ' ', ljac$hh)
ljac$date.time <- strptime(ljac$date.time, format = '%Y-%m-%d %H', tz = 'GMT')
ljac$date.time <- as.POSIXlt(ljac$date.time, tz = 'PST')
ljac.hourly <- aggregate(ljac[,ljac.col.select], by = list(as.character(ljac$date.time)), FUN = mean, na.rm = T)

### NOAA - ljpc1 ###

## https://www.ndbc.noaa.gov/station_history.php?station=ljpc1

ljpc1.col.select <- c('WVHT','DPD')
ljpc1 <- read.table('NOAA_ljpc1/combined_ljpc.txt', col.names = noaa.columns, na.strings = noaa.na)
ljpc1$date.time <- paste0(ljpc1$YY, '-', ljpc1$MM, '-', ljpc1$DD, ' ', ljpc1$hh)
ljpc1$date.time <- strptime(ljpc1$date.time, format = '%Y-%m-%d %H', tz = 'GMT')
ljpc1$date.time <- as.POSIXlt(ljpc1$date.time, tz = 'PST')
ljpc1.hourly <- aggregate(ljpc1[,ljpc1.col.select], by = list(as.character(ljpc1$date.time)), FUN = mean, na.rm = T)

#### combine datasets ####

library(zoo)

mims.zoo <- zoo(mims.hourly, order.by = mims.hourly$Group.1)
miniDot.zoo <- zoo(miniDot.hourly, order.by = miniDot.hourly$Group.1)
sccoos.zoo <- zoo(sccoos.hourly, order.by = sccoos.hourly$Group.1)
ljac.zoo <- zoo(ljac.hourly, order.by = ljac.hourly$Group.1)
ljpc1.zoo <- zoo(ljpc1.hourly, order.by = ljpc1.hourly$Group.1)

saveRDS(mims.zoo, file = "2023-02-03_mims_zoo.rds")
saveRDS(miniDot.zoo, file = "2023-02-03_miniDot_zoo.rds")
saveRDS(sccoos.zoo, file = "2023-02-03_sccoos_zoo.rds")
saveRDS(ljac.zoo, file = "2023-02-03_ljac_zoo.rds")
saveRDS(ljpc1.zoo, file = "2023-02-03_ljpc1_zoo.rds")


combined.zoo <- merge(mims.zoo, miniDot.zoo)
combined.zoo <- merge(combined.zoo, sccoos.zoo)
combined.zoo <- merge(combined.zoo, ljac.zoo)
combined.zoo <- merge(combined.zoo, ljpc1.zoo)

combined.df <- as.data.frame(combined.zoo)
combined.df <- as.data.frame(lapply(combined.df, as.numeric))
combined.df$dates <- strptime(index(combined.zoo), format = '%Y-%m-%d %H')

## The aou calculation is done here, because you either need the O2_sat value
## that comes with the MIMS data or salinity which comes from SCCOOS

combined.df$O2_sat <- O2sat(t = combined.df$Temperature)
combined.df$aou <- -1 * (combined.df$O2_sat - combined.df$Dissolved.Oxygen)
combined.df$delta <- combined.df$o2_bio - combined.df$aou

### some quality checks, could probably do this earlier ###

hist(combined.df$ATMP)
hist(combined.df$PRES)
hist(combined.df$DPD)
hist(combined.df$WVHT)
hist(combined.df$GST)
hist(combined.df$WSPD)
hist(combined.df$pressure)
hist(combined.df$salinity, breaks = 100)
hist(combined.df$temperature.sccoos.zoo)
hist(combined.df$Dissolved.Oxygen, breaks = 1000)

plot(combined.df$dates, combined.df$Dissolved.Oxygen, type = 'l')

plot(combined.df$o2_bio ~ combined.df$aou)
plot(combined.df$dates, combined.df$aou)

plot(combined.df$dates, combined.df$aou)

points(combined.df$dates, combined.df$o2_bio,
       col = 'red')

## salinity has many problematic values, exclude bad time points
plot(combined.df$salinity~as.Date(combined.df$dates))

combined.df <- combined.df[which(combined.df$salinity > 33.3 & combined.df$salinity < 33.9),]

plot(combined.df$salinity~as.Date(combined.df$dates))

saveRDS(combined.df, file = "2023-03-08_combined_env_data_hourly.rds")


## aggregate by day

# library(lubridate)
# library(tidyverse)
# combined.df$Date <- parse_date_time(paste(day(combined.df$dates), month(combined.df$dates), year(combined.df$dates), sep = "-"), orders = "dmY")
# my.colnames <- colnames(combined.df)[c(1:27, 29:30)]
# my.cols <- c(1:27, 29:30)
# combined.df.daily <- combined.df %>% group_by(Date) %>% summarize_at(vars(my.colnames), mean)
# saveRDS(combined.df.daily, file = "2022-10-10_combined_env_data_daily.rds")


## nice plots of available data

plot(combined.df$dates, combined.df$aou,
     type = 'l',
     ylab = expression("[O"[2]*"]"[bio]*"|AOU  ["*mu*"M]"),
     xlab = 'Date')

points(combined.df$dates, combined.df$o2_bio, type = 'l', col = 'blue')

legend('topleft',
       legend = c('AOU', expression("[O"[2]*"]"[bio])),
       lty = 1,
       col = c('black', 'blue'))

# # just AOU
# plot(combined.df$dates, combined.df$aou,
#      type = 'l',
#      ylab = expression("AOU  ["*mu*"M]"),
#      xlab = 'Date')
# 
# legend('topleft',
#        legend = c('AOU'),
#        lty = 1,
#        col = c('black'))
# 
# # just [o2]bio
# plot(combined.df$dates, combined.df$o2_bio,
#      type = 'l',
#      col = "blue",
#      ylab = expression("[O"[2]*"]"[bio]*"  ["*mu*"M]"),
#      xlab = 'Date')
# 
# legend('topleft',
#        legend = c(expression("[O"[2]*"]"[bio])),
#        lty = 1,
#        col = c('blue'))


#### train and validate boosted regression tree ####

predictors <- c('aou',
                #                'delta',
                "temperature.sccoos.zoo",
                'salinity',
                "pressure",
                "WSPD",
                "GST",
                "WVHT",
                "DPD",
                #                "ATMP",
                #"PRES",
                "o2_bio")

accessory <- c('dates', 'Dissolved.Oxygen', 'O2')

## check to see where O2 and O2 agree

plot(combined.df$O2[combined.df$dates < strptime('2021-3-26', format = '%Y-%m-%d')],
     combined.df$Dissolved.Oxygen[combined.df$dates < strptime('2021-3-26', format = '%Y-%m-%d')],
     log = 'x')

plot(combined.df$O2,
     combined.df$Dissolved.Oxygen,
     log = 'x')

## data selection for modeling, constrain by vacuum pressure (MIMS)

combined.df.select <- combined.df[which(combined.df$O2 > 1.5e-9), c(predictors, accessory)]

## not sure why there are missing values, but seem to be the case!

combined.df.select <- na.omit(combined.df.select)

saveRDS(combined.df.select, file = "2023-03-08_combined.df.select.rds")

## date range for validation

date1 <- strptime('2021-7-1', format = '%Y-%m-%d')
date2 <- strptime('2021-7-15', format = '%Y-%m-%d')

combined.df.train <- combined.df.select

combined.df.test <- combined.df.select[which(combined.df.select$dates > date1 &
                                               combined.df.select$dates < date2),]

combined.df.train <- combined.df.train[which(!combined.df.train$dates %in% combined.df.test$dates),]

combined.df.test.dates <- combined.df.test$dates
combined.df.train.dates <- combined.df.train$dates

combined.df.train <- combined.df.train[,predictors]
combined.df.test <- combined.df.test[,predictors]

plot(combined.df.select$dates, combined.df.select$o2_bio, col = 'blue', type = 'l',
     ylab = expression("[O"[2]*"]"[bio]*"|AOU  ["*mu*"M]"),
     xlab = "Date")
points(combined.df.select$dates, combined.df.select$aou, type = 'l')
points(combined.df.test.dates, combined.df.test$o2_bio, col = 'red', type = 'l')

legend('topleft',
       legend = c('AOU', expression("[O"[2]*"]"[bio]), 'Testing'),
       col = c('black', 'blue', 'red'),
       lty = 1)

## train partial model

library(gbm)

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

plot(combined.gbm$train.error)

combined.gbm.summary <- summary(combined.gbm) # assess variable significance

# 53.6+29.7+6.5

plot(combined.gbm$fit ~ combined.df.train$o2_bio,
     ylab = 'Predicted',
     xlab = 'Observed')

abline(0, 1, col = 'red')

summary(lm(combined.gbm$fit ~ combined.df.train$o2_bio))

plot((combined.gbm$fit - combined.df.train$o2_bio),
     type = 'h')

#### apply model to validation data ####

test.aou <- predict(combined.gbm, combined.df.test)

test.aou.o2bio.delta <- test.aou - combined.df.test$o2_bio #difference between predicted and o2bio
aou.o2bio.delta <- combined.df.test$aou - combined.df.test$o2_bio #difference between aou and o2bio
aou.o2bio.rmse <- mean(test.aou.o2bio.delta**2)

plot(combined.df.test.dates,
     aou.o2bio.delta,
     type = 'h')

points(combined.df.test.dates,
       test.aou.o2bio.delta,
     type = 'h',
     col = 'red')

hist(test.aou.o2bio.delta, breaks = 100, col = 'black') # corrected
hist(aou.o2bio.delta, breaks = 100, col = 'red', add = T) # not corrected

plot(combined.df.test.dates,
     test.aou,
     type = 'l',
     ylab = expression("[O"[2]*"]"[bio]*"|AOU  ["*mu*"M]"),
     ylim = c(-100, 100),
     col = 'orange',
     xlab = 'Date')

points(combined.df.test.dates,
       (combined.df.test$aou),
       type = 'l',
       col = 'black')

points(combined.df.test.dates,
       combined.df.test$o2_bio,
       type = 'l',
       col = 'blue')

legend('topleft',
       legend = c('AOU', 'AOU corrected', expression("[O"[2]*"]"[bio])),
       col = c('black', 'orange', 'blue'),
       lty = 1)






plot(test.aou ~ combined.df.test$o2_bio,
     col = 'blue')
points(combined.df.test$aou ~ combined.df.test$o2_bio,
       col = 'black')
abline(1,1)
summary(lm(test.aou ~ combined.df.test$o2_bio))
summary(lm(combined.df.test$aou ~ combined.df.test$o2_bio))

#### parameter hypertuning ####

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
    hyper.grid$RMSE[i] <- mean(temp.delta**2)
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

# saveRDS(final.gbm, file = "2022-10-03_final_gbm.rds")
final.gbm <- readRDS("2022-10-03_final_gbm.rds")

## apply final model to test data

final.test.aou <- predict(final.gbm, combined.df.test)

plot(combined.df.test.dates,
     final.test.aou,
     type = 'l',
     ylab = expression("[O"[2]*"]"[bio]*"|AOU  ["*mu*"M]"),
     ylim = c(-100, 100),
     col = 'orange',
     xlab = 'Date')

points(combined.df.test.dates,
       (combined.df.test$aou),
       type = 'l',
       col = 'black')

points(combined.df.test.dates,
       combined.df.test$o2_bio,
       type = 'l',
       col = 'blue')

legend('topleft',
       legend = c('AOU', 'AOU corrected', expression("[O"[2]*"]"[bio])),
       col = c('black', 'orange', 'blue'),
       lty = 1)

#### full model ####

full.gbm <- gbm(o2_bio ~ .,
                    data = na.omit(combined.df[,predictors]),
                    interaction.depth = hyper.grid$interaction.depth[which.min(hyper.grid$RMSE)],
                    distribution = 'gaussian',
                    n.trees       = hyper.grid$n.trees[which.min(hyper.grid$RMSE)],
                    bag.fraction  = hyper.grid$bag.fraction[which.min(hyper.grid$RMSE)],
                    shrinkage   = hyper.grid$shrinkage[which.min(hyper.grid$RMSE)],
                    cv.folds = 2)

save(list = c('combined.gbm', 'full.gbm'), file = '20220930_model.Rdata')

#### apply model to miniDot timeseries ####

load(file = "20220930_model.Rdata")

full.predictors <- na.omit(combined.df[,c(full.gbm$var.names, 'dates')])

full.aou.corrected <- predict(full.gbm, full.predictors)

plot(full.predictors$dates, full.aou.corrected, type = 'l')
points(full.predictors$dates, full.predictors$aou, type = 'l', col = 'red')
points(combined.df.select$dates, combined.df.select$o2_bio, type = 'l', col = 'blue')

## code for smoothing functions needs to be corrected from full.df (no longer exists)
## to combined.df.select

smooth.i <- 25 
smooth.ends <- floor((smooth.i) / 2)

full.aou.corrected.smooth <- rollmean(full.aou.corrected, smooth.i)
full.aou.smooth <- rollmean(full.predictors$aou, smooth.i)
full.o2bio.smooth <- rollmean(combined.df$o2_bio, smooth.i)

plot(combined.df$dates[(smooth.ends + 1):(dim(combined.df)[1] - smooth.ends)], full.o2bio.smooth, type = 'l', col = 'blue',
     ylab = expression("AOU  ["*mu*"M]"),
     xlab = "Date")
points(full.predictors$dates[(smooth.ends + 1):(dim(full.predictors)[1] - smooth.ends)], full.aou.smooth, type = 'l', col = 'black')
points(full.predictors$dates[(smooth.ends + 1):(dim(full.predictors)[1] - smooth.ends)], full.aou.corrected.smooth, type = 'l', col = "orange")

legend('topleft',
       legend = c('AOU', 'AOU_corrected', '[O2]bio'),
       col = c('black', 'orange', 'blue'),
       lty = 1)

full.aou.df <- as.data.frame(cbind(full.aou.corrected, full.predictors$aou))
colnames(full.aou.df) <- c('aou.corrected', 'aou')
full.aou.df$dates <- full.predictors$dates

full.aou.df <- merge(full.aou.df, combined.df, by = c("dates", "aou"))

write.csv(full.aou.df, 'corrected_aou.csv')

#### cartoon model ####

## This is a toy decision tree that is useful for talking to an audience about
## gradient boosted trees, etc.

dummy.predictors <- c("o2_bio",
"temperature.sccoos.zoo",
                'salinity',
                "pressure",
                "WSPD",
                "GST",
                "WVHT",
                "DPD")
library(tree)

temp <- tree(o2_bio ~ .,
             data = na.omit(combined.df[,predictors]))
plot(temp)
text(temp)


#### export data ####

saveRDS(full.aou.df, "full_aou_df.rds")
saveRDS(full.predictors, "full_predictors_df.rds")
saveRDS(mims.hourly, "mims_hourly.rds")


#### to do items ####

## What drive the difference between AOU and AOU_corrected?  This can be explored in a number of ways, ideally through an ordination (e.g. PCA)
## of the predictors following proper normalization.  Then, look to see if AOU-AOU_corrected correlates to any of the principal components.