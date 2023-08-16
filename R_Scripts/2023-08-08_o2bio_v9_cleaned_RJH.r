## A more correct approach would be to derive a "true" O2 value for the MIMS
## based on existing calibrations, then base the model only on timepoints where
## the mims and miniDOT O2 values agree
getwd()
setwd('C://Users/haler/Documents/PhD-Bowman/O2-Ar_time_series/R_Data/')
set.seed(1234)


# ---- library ----

library(tidyverse)
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




# ---- get data ----

## Get O2bio calculated in read_lvm.py ##

mims <- read.csv('Current_Model_Inputs/MIMS_o2bio/o2bio.csv') # this file is produced by read_lvm.py, refer to that script for details
mims2 <- read.csv('Current_Model_Inputs/MIMS_o2bio/o2bio_vol2.csv')
mims <- rbind(mims, mims2) # combine datasets
mims$date_time <- strptime(mims$date_time, format = '%Y-%m-%d %H')
mims1 <- mims[mims$N2.Ar < 40 & mims$N2.Ar > 30 & mims$date_time < strptime('2021-3-26', format = '%Y-%m-%d'),] # manual QC
mims2 <- mims[mims$N2.Ar < 20 & mims$N2.Ar > 9 & mims$date_time >= strptime('2021-3-26', format = '%Y-%m-%d'),]
mims <- rbind(mims1, mims2)
mims.date <- mims$date_time

## aggregate hourly
mims.hourly <- mims %>% group_by(date_time) %>% 
  summarize(temperature = mean(temperature), O2_sat = mean(O2_sat), Ar_sat = mean(Ar_sat), O2.Ar_sat = mean(O2.Ar_sat), 
            O2 = mean(O2), O2.Ar = mean(O2.Ar), N2.Ar = mean(N2.Ar), O2_CF = mean(O2_CF), o2_bio = mean(o2_bio))

## quick data visualization
plot(mims.date, mims$o2_bio)
plot(mims.date, mims$N2.Ar)



### Get miniDOT data ###

miniDot <- read.csv("Current_Model_Inputs/SIOPier_miniDOT_20180809_20230131_JBowman.csv", header = F, na.strings = "NA", skip = 2)
miniDot <- miniDot[,-1]
colnames(miniDot) <- c("UTC_Date_._Time", "Pacific.Standard.Time", "Temperature", "Dissolved.Oxygen", "Dissolved.Oxygen.Saturation", "Sensor")

miniDot$Pacific.Standard.Time <- strptime(miniDot$Pacific.Standard.Time, format = '%Y-%m-%d %H')

## MiniDot reports in mg/L, need uMol
miniDot$Dissolved.Oxygen <- (miniDot$Dissolved.Oxygen / (15.999 * 2 * 1000)) * 10 ** 6

miniDot.col.select <- c("Dissolved.Oxygen", "Temperature")

## aggregate hourly ##
miniDot.hourly <- miniDot %>% group_by(Pacific.Standard.Time) %>% 
  summarize(Temperature = mean(Temperature), Dissolved.Oxygen = mean(Dissolved.Oxygen), Dissolved.Oxygen.Saturation = mean(Dissolved.Oxygen.Saturation))


### Get sccoos data ###

sccoos <- read.csv('https://erddap.sccoos.org/erddap/tabledap/autoss.csv?station%2Ctime%2Ctemperature%2Cconductivity%2Cpressure%2Cchlorophyll%2Csalinity&station=%22scripps_pier%22')
sccoos <- sccoos[-1,]
sccoos.columns <- c('station', 'time', 'temperature', 'conductivity','pressure','chlorophyll','salinity')
sccoos$station <- NULL
sccoos$time <- strptime(sccoos$time, format = '%Y-%m-%dT%H', tz = 'GMT')
sccoos$time <- as.POSIXlt(sccoos$time, tz = 'PST')
sccoos$temperature <- as.numeric(sccoos$temperature)
sccoos$pressure <- as.numeric(sccoos$pressure)
sccoos$chlorophyll <- as.numeric(sccoos$chlorophyll)
sccoos$conductivity <- as.numeric(sccoos$conductivity)

## aggregate hourly ##
sccoos.hourly <- sccoos %>% group_by(time) %>% 
  summarize(temperature.sccoos = mean(temperature), conductivity = mean(conductivity), pressure = mean(pressure), chlorophyll = mean(chlorophyll), salinity = mean(salinity))


## can aggregate daily if desired ##

# sccoos.daily <- sccoos
# sccoos.daily$sample.date <- parse_date_time(paste(year(sccoos.daily$time), month(sccoos.daily$time), day(sccoos.daily$time), sep = "-"), orders = "Ymd")
# sccoos.daily <- sccoos.daily %>% group_by(sample.date) %>%
#   summarize(temperature.sccoos = mean(temperature), conductivity = mean(conductivity), pressure = mean(pressure), chlorophyll = mean(chlorophyll), salinity = mean(salinity))
# saveRDS(sccoos.daily, file = "2023-04-12_sccoos_env_data_daily_mean.rds")


### Get NOAA LJAC1 buoy data ###
## Data downloaded from: https://www.ndbc.noaa.gov/station_history.php?station=ljac1

noaa.columns <- c('YY','MM','DD','hh','mm','WDIR','WSPD','GST','WVHT','DPD','APD','MWD','PRES','ATMP','WTMP','DEWP','VIS','TIDE')
noaa.na <- c('99.00', '999', '999.0', '9999.0', '99.0')

ljac.col.select <- c('WSPD','GST','PRES','ATMP','WTMP')

# files all saved in same folder
file.names <- list.files(path = "Current_Model_Inputs/NOAA_LJAC1_RJH/")
for (i in seq_along(file.names)) {
  assign(substr(file.names[i], start = 1, stop = nchar(file.names[i])-4), read.table(paste("Current_Model_Inputs/NOAA_LJAC1_RJH/", file.names[i], sep = ""), col.names = noaa.columns, na.strings = noaa.na))
}

ljac <- rbind(ljac1h2018, ljac1h2019, ljac1h2020, ljac1h2021, ljac1h2022, ljac1jan2023, ljac1feb2023, ljac1mar2023, ljac1apr2023, ljac1may2023, ljac1jun2023)
ljac$date.time <- paste0(ljac$YY, '-', ljac$MM, '-', ljac$DD, ' ', ljac$hh)
ljac$date.time <- strptime(ljac$date.time, format = '%Y-%m-%d %H', tz = 'GMT')
ljac$date.time <- as.POSIXlt(ljac$date.time, tz = 'PST')

## aggregate hourly
ljac.hourly <- ljac %>% group_by(date.time) %>% 
  summarize(WSPD = mean(WSPD), GST = mean(GST), PRES = mean(PRES), ATMP = mean(ATMP), WTMP = mean(WTMP))


### Get NOAA LJPC1 buoy data ###
## Data downloaded from: https://www.ndbc.noaa.gov/station_history.php?station=ljpc1

noaa.columns <- c('YY','MM','DD','hh','mm','WDIR','WSPD','GST','WVHT','DPD','APD','MWD','PRES','ATMP','WTMP','DEWP','VIS','TIDE')
noaa.na <- c('99.00', '999', '999.0', '9999.0', '99.0')

ljpc1.col.select <- c('WVHT','DPD')

# files all saved in same folder
file.names <- list.files(path = "Current_Model_Inputs/NOAA_LJPC1_RJH/")
for (i in seq_along(file.names)) {
  assign(substr(file.names[i], start = 1, stop = nchar(file.names[i])-4), read.table(paste("Current_Model_Inputs/NOAA_LJPC1_RJH/", file.names[i], sep = ""), col.names = noaa.columns, na.strings = noaa.na))
}

ljpc1 <- rbind(ljpc1h2018, ljpc1h2019, ljpc1h2020, ljpc1h2021, ljpc1h2022, ljpc1jan2023, ljpc1feb2023, ljpc1mar2023, ljpc1apr2023, ljpc1may2023, ljpc1jun2023)
ljpc1$date.time <- paste0(ljpc1$YY, '-', ljpc1$MM, '-', ljpc1$DD, ' ', ljpc1$hh)
ljpc1$date.time <- strptime(ljpc1$date.time, format = '%Y-%m-%d %H', tz = 'GMT')
ljpc1$date.time <- as.POSIXlt(ljpc1$date.time, tz = 'PST')

## aggregate hourly
ljpc1.hourly <- ljpc1 %>% group_by(date.time) %>% 
  summarize(WVHT = mean(WVHT), DPD = mean(DPD))



# ---- combine datasets ----

colnames(mims.hourly)[1] <- "Date.Time"
colnames(miniDot.hourly)[1] <- "Date.Time"
colnames(sccoos.hourly)[1] <- "Date.Time"
colnames(ljac.hourly)[1] <- "Date.Time"
colnames(ljpc1.hourly)[1] <- "Date.Time"

mims.hourly <- data.frame(mims.hourly)
miniDot.hourly <- data.frame(miniDot.hourly)
sccoos.hourly <- data.frame(sccoos.hourly)
ljac.hourly <- data.frame(ljac.hourly)
ljpc1.hourly <- data.frame(ljpc1.hourly)

mims.hourly$Date.Time <- as.character(mims.hourly$Date.Time)
miniDot.hourly$Date.Time <- as.character(miniDot.hourly$Date.Time)
sccoos.hourly$Date.Time <- as.character(sccoos.hourly$Date.Time)
ljac.hourly$Date.Time <- as.character(ljac.hourly$Date.Time)
ljpc1.hourly$Date.Time <- as.character(ljpc1.hourly$Date.Time)

combined.df <- merge(mims.hourly, miniDot.hourly, by = "Date.Time", all = T, sort = T)
combined.df <- merge(combined.df, sccoos.hourly, by = "Date.Time", all.x = T, sort = T)
combined.df <- merge(combined.df, ljac.hourly, by = "Date.Time", all.x = T, sort = T)
combined.df <- merge(combined.df, ljpc1.hourly, by = "Date.Time", all.x = T, sort = T)

## turn date into PociXct format using parse_date_time() function
combined.df$Date.Time <- parse_date_time(combined.df$Date.Time, orders = "Ymd HMS")


# ---- calculate AOP, O2_sat, and delta ----

## The AOP calculation is done here, because you either need the O2_sat value
## that comes with the MIMS data or salinity which comes from SCCOOS
## Note: AOP = -AOU

combined.df$O2_sat <- O2sat(t = combined.df$Temperature, s = combined.df$salinity) # uses function defined at beginning
combined.df$aop <- -1 * (combined.df$O2_sat - combined.df$Dissolved.Oxygen)
combined.df$delta <- combined.df$o2_bio - combined.df$aop



# ---- QC data ----

hist(combined.df$ATMP)
hist(combined.df$PRES)
hist(combined.df$DPD)
hist(combined.df$WVHT)
hist(combined.df$GST)
hist(combined.df$WSPD)
hist(combined.df$pressure)
hist(combined.df$salinity, breaks = 100)
hist(combined.df$temperature.sccoos)
hist(combined.df$Dissolved.Oxygen, breaks = 1000)

plot(combined.df$Date.Time, combined.df$Dissolved.Oxygen, type = 'l')

plot(combined.df$o2_bio ~ combined.df$aop)
plot(combined.df$Date.Time, combined.df$aop)

## compare MIMS o2bio time series and miniDOT AOP time series
plot(combined.df$Date.Time, combined.df$aop)
points(combined.df$Date.Time, combined.df$o2_bio,
       col = 'red')


## salinity has many problematic values, exclude bad time points
plot(combined.df$salinity~as.Date(combined.df$Date.Time))

combined.df <- combined.df[which(combined.df$salinity > 33.3 & combined.df$salinity < 33.9),]

plot(combined.df$salinity~as.Date(combined.df$Date.Time))

saveRDS(combined.df, file = "2023-08-08_combined_env_data_hourly.rds")


# ---- aggregate by day ----

combined.df$Date <- parse_date_time(paste(day(combined.df$Date.Time), month(combined.df$Date.Time), year(combined.df$Date.Time), sep = "-"), orders = "dmY")
my.colnames <- colnames(combined.df)[c(2:27)]
combined.df.daily <- combined.df %>% group_by(Date) %>% summarize_at(vars(my.colnames), mean)
saveRDS(combined.df.daily, file = "2023-08-08_combined_env_data_daily.rds")


# ---- nice plots of available data -----

plot(combined.df$Date, combined.df$aop,
     type = 'l',
     ylab = expression("[O"[2]*"]"[bio]*"|AOP  ["*mu*"M]"),
     xlab = 'Date')
points(combined.df$Date, combined.df$o2_bio, type = 'l', col = 'blue')
legend('topleft',
       legend = c('AOU', expression("[O"[2]*"]"[bio])),
       lty = 1,
       col = c('black', 'blue'))


ggplot(data = combined.df) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  geom_line(aes(x = Date.Time, y = aop), color = col1.aou, lwd = 1, alpha = 0.7) +
  geom_line(aes(x = Date.Time, y = o2_bio), color = col2.o2bio, lwd = 1, alpha = 0.7) +
  labs(x = "Date", y = expression("Oxygen Anomaly  [ "*mu*"M]")) +
  theme_bw() +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  theme(panel.grid = element_blank()) +
  ylim(c(-250,250)) +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12))



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


## check to see where O2 and O2 agree ##

plot(combined.df$O2[combined.df$Date.Time < strptime('2021-3-26', format = '%Y-%m-%d')],
     combined.df$Dissolved.Oxygen[combined.df$Date.Time < strptime('2021-3-26', format = '%Y-%m-%d')],
     log = 'x')

plot(combined.df$O2,
     combined.df$Dissolved.Oxygen,
     log = 'x')


## data selection for modeling, constrain by vacuum pressure (MIMS) ##
combined.df.select <- combined.df[which(combined.df$O2 >= 1.5e-9), c(predictors, accessory)]

## not sure why there are missing values, but seem to be the case!
combined.df.select <- na.omit(combined.df.select)

saveRDS(combined.df.select, file = "2023-08-08_combined.df.select.rds")


sqrt(mean((combined.df.select$o2_bio - combined.df.select$aop)^2))




## date range for validation

date1 <- strptime('2021-7-1', format = '%Y-%m-%d')
date2 <- strptime('2021-7-15', format = '%Y-%m-%d')

combined.df.train <- combined.df.select

combined.df.test <- combined.df.select[which(combined.df.select$Date.Time > date1 &
                                               combined.df.select$Date.Time < date2),]

combined.df.train <- combined.df.train[which(!combined.df.train$Date.Time %in% combined.df.test$Date.Time),]

combined.df.test.Date.Time <- combined.df.test$Date.Time
combined.df.train.Date.Time <- combined.df.train$Date.Time

combined.df.train <- combined.df.train[,predictors]
combined.df.test <- combined.df.test[,predictors]


ggplot() +
  geom_hline(aes(yintercept = 0), alpha = 0.2) +
  geom_line(aes(x = combined.df.select$Date.Time[-c(1:2)], y = combined.df.select$o2_bio[-c(1:2)]), color = col2.o2bio, lwd = 1, alpha = 0.7) +
  geom_line(aes(x = combined.df.select$Date.Time[-c(1:2)], y = combined.df.select$aop[-c(1:2)]), color = col1.aou, lwd = 1, alpha = 0.7) +
  # geom_line(aes(x = combined.df.test.Date.Time, y = combined.df.test$o2_bio), color = col3, lwd = 1, alpha = 1) +
  labs(x = "Date", y = expression("[O"[2]*"]"[bio]*" | AOP  ["*mu*"M]")) +
  ylim(c(-250,250)) +
  theme_bw() +
  scale_x_datetime(date_breaks = "3 months", date_labels = "%b %Y") +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12))
  



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

plot(combined.gbm$train.error)

## assess variable significance ##
combined.gbm.summary <- summary(combined.gbm) 

combined.gbm.summary$var <- factor(combined.gbm.summary$var, levels = rev(rownames(combined.gbm.summary)))
ggplot(data = combined.gbm.summary) +
  geom_col(aes(x = var, y = rel.inf), fill = c(col5.other, col5.other, col5.other, col1.aou, col1.aou, col1.aou, col1.aou, col1.aou)) +
  # scale_fill_manual(values = ) +
  labs(x = "GBM Predictor", y = "Relative Influence") +
  scale_x_discrete(labels = rev(c("AOP", "Water Temperature", "Salinity", "Wave Height", "Dominant Period", "Gust Speed", "Pressure", "Wind Speed"))) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12))

sum(combined.gbm.summary$rel.inf[1:3])

## assess model fit ##
plot(combined.gbm$fit ~ combined.df.train$o2_bio,
     ylab = 'Predicted',
     xlab = 'Observed')

abline(0, 1, col = 'red')

summary(lm(combined.gbm$fit ~ combined.df.train$o2_bio))

plot((combined.gbm$fit - combined.df.train$o2_bio),
     type = 'h')



# ---- apply model to validation data ----

aop.cor <- predict(combined.gbm, combined.df.test)

aop.cor.delta <- aop.cor - combined.df.test$o2_bio #difference between predicted and o2bio
aop.o2bio.delta <- combined.df.test$aop - combined.df.test$o2_bio #difference between uncorrected aop and o2bio

aop.cor.rmse <- sqrt(mean(aop.cor.delta^2))

plot(combined.df.test.Date.Time,
     aou.o2bio.delta,
     type = 'h')

points(combined.df.test.Date.Time,
       aop.cor.delta,
     type = 'h',
     col = 'red')

hist(aop.cor.delta, breaks = 100, col = 'black') # corrected
hist(aop.o2bio.delta, breaks = 100, col = 'red', add = T) # not corrected

plot(combined.df.test.Date.Time,
     test.aop,
     type = 'l',
     ylab = expression("[O"[2]*"]"[bio]*"|AOU  ["*mu*"M]"),
     ylim = c(-100, 100),
     col = 'orange',
     xlab = 'Date')

points(combined.df.test.Date.Time,
       (combined.df.test$aop),
       type = 'l',
       col = 'black')

points(combined.df.test.Date.Time,
       combined.df.test$o2_bio,
       type = 'l',
       col = 'blue')

legend('topleft',
       legend = c('AOU', 'AOU corrected', expression("[O"[2]*"]"[bio])),
       col = c('black', 'orange', 'blue'),
       lty = 1)


plot(test.aop ~ combined.df.test$o2_bio,
     col = 'blue')
points(combined.df.test$aop ~ combined.df.test$o2_bio,
       col = 'black')
abline(1,1)
summary(lm(test.aop ~ combined.df.test$o2_bio))
summary(lm(combined.df.test$aop ~ combined.df.test$o2_bio))



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

# saveRDS(final.gbm, file = "2023-08-08_final_gbm.rds")
final.gbm <- readRDS("2023-08-08_final_gbm.rds")

## apply final model to test data

aop.cor.final <- predict(final.gbm, combined.df.test)

plot(combined.df.test.Date.Time,
     aop.cor.final,
     type = 'l',
     ylab = expression("[O"[2]*"]"[bio]*"|AOP  ["*mu*"M]"),
     ylim = c(-100, 100),
     col = 'orange',
     xlab = 'Date')

points(combined.df.test.Date.Time,
       (combined.df.test$aop),
       type = 'l',
       col = 'black')

points(combined.df.test.Date.Time,
       combined.df.test$o2_bio,
       type = 'l',
       col = 'blue')

legend('topleft',
       legend = c('AOP', 'AOP corrected', expression("[O"[2]*"]"[bio])),
       col = c('black', 'orange', 'blue'),
       lty = 1)

ggplot() +
  geom_hline(aes(yintercept = 0), alpha = 0.2) +
  geom_line(aes(x = combined.df.test.Date.Time, y = combined.df.test$o2_bio), color = col2.o2bio, lwd = 1, alpha = 0.7) +
  geom_line(aes(x = combined.df.test.Date.Time, y = combined.df.test$aop), color = col1.aou, lwd = 1, alpha = 0.7) +
  geom_line(aes(x = combined.df.test.Date.Time, y = aop.cor.final), color = col3.aoucor, lwd = 1, alpha = 0.7) +
  labs(x = "Date", y = expression("Oxygen Anomaly  [ "*mu*"M]")) +
  ylim(c(-120,70)) +
  theme_bw() +
  scale_x_datetime(date_breaks = "3 days", date_labels = "%b %d") +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12)) 


# ---- calculate model performance ----

my.rmse <- sqrt(mean((combined.df.test$o2_bio - aop.cor.final)^2))
my.rmse.old <- sqrt(mean((combined.df.test$o2_bio - combined.df.test$aop)^2))

# cor(combined.df.test$o2_bio, final.test.aou)
# cor(combined.df.test$o2_bio, combined.df.test$aou)

model <- lm(aop.cor.final~combined.df.test$o2_bio)
model.old <- lm(combined.df.test$aop~combined.df.test$o2_bio)

summary(model)
summary(model.old)

# wald.test(vcov(model), b = coef(model), Terms = 1, H0 = 1)
# wald.test(vcov(model.old), b = coef(model), Terms = 1, H0 = 1)

# t.stat <- (coef(summary(model))[2,1] - 1)/(coef(summary(model))[2,1])

t.test(x = (aop.cor.final-combined.df.test$o2_bio), mu = 0)
t.test(x = (combined.df.test$aop-combined.df.test$o2_bio), mu = 0)

hist(aop.cor.final-combined.df.test$o2_bio)
hist(combined.df.test$aop-combined.df.test$o2_bio)


sum(abs(aop.cor.final-combined.df.test$o2_bio))
sum(abs(combined.df.test$aop-combined.df.test$o2_bio))



ggplot() +
  geom_point(aes(x = combined.df.test$o2_bio, y = combined.df.test$aop), color = col1.aou, size = 2) +
  geom_point(aes(combined.df.test$o2_bio, y = aop.cor.final), color = col3.aoucor, size = 2) +
  geom_abline(slope = 1, intercept = 0, color = col2.o2bio, size = 2) +
  scale_y_continuous(breaks = c(-75, -50, -25, 0, 25)) +
  labs(x = "[O2]bio", y = "AOP | Corrected AOP") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12))


# ---- full model ----

full.gbm <- gbm(o2_bio ~ .,
                    data = na.omit(combined.df[,predictors]),
                    interaction.depth = hyper.grid$interaction.depth[which.min(hyper.grid$RMSE)],
                    distribution = 'gaussian',
                    n.trees       = hyper.grid$n.trees[which.min(hyper.grid$RMSE)],
                    bag.fraction  = hyper.grid$bag.fraction[which.min(hyper.grid$RMSE)],
                    shrinkage   = hyper.grid$shrinkage[which.min(hyper.grid$RMSE)],
                    cv.folds = 2)

save(list = c('combined.gbm', 'full.gbm'), file = '20230808_model.Rdata')



# ---- apply model to miniDot timeseries ----

combined.df <- readRDS("2023-08-08_combined_env_data_hourly.rds")
load(file = "20230808_model.Rdata")

full.predictors <- na.omit(combined.df[,c(full.gbm$var.names, 'Date.Time')])

full.predictors$aop.corrected <- predict(full.gbm, full.predictors)

saveRDS(full.predictors, "2023-08-08_aop_cor_df.rds")

ggplot() +
  geom_hline(aes(yintercept = 0), alpha = 0.5) +
  geom_line(aes(x = full.predictors$Date.Time, y = full.predictors$aop), color = col1.aou, lwd = 1, alpha = 0.7) +
  geom_line(aes(x = combined.df.select$Date.Time, y = combined.df.select$o2_bio), color = col2.o2bio, lwd = 1, alpha = 0.7) +
  geom_line(aes(x = full.predictors$Date.Time, y = full.aop.corrected), color = col3.aoucor, lwd = 1, alpha = 0.7) +
  labs(x = "Date", y = expression("Oxygen Anomaly  [ "*mu*"M]")) +
  theme_bw() +
  ylim(c(-250, 250)) +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12)) 


# ---- smoothing ----

date.start <- full.predictors$Date.Time[1] - weeks(1)
date.end <- full.predictors$Date.Time[nrow(full.predictors)] + weeks(1)
test <- seq(from = date.start, to = date.end, by = "hours")

full.predictors.with.aop.cor <- full.predictors
full.predictors.with.aop.cor$aou_corrected <- full.aop.corrected

my.rolling.avgs.df <- full.predictors.with.aop.cor

# clunky but it works
for(r in 1:nrow(full.predictors.with.aop.cor)){
  
  temp.df <- full.predictors.with.aop.cor[which(full.predictors.with.aop.cor$Date.Time >= (full.predictors.with.aop.cor$Date.Time[r] - weeks(1)) & full.predictors.with.aop.cor$Date.Time <= (full.predictors.with.aop.cor$Date.Time[r] + weeks(1))),]
  
  my.rolling.avgs.df[r,-9] <- colMeans(temp.df[-9])
  
}



sccoos.dates <- readRDS("../16S_sccoos/2023-04-28_sccoos_dates.rds")

ggplot() +
  geom_line(aes(x = full.predictors$Date.Time, y = full.aop.corrected), color = col3.aoucor, lwd = 1, alpha = 0.7) +
  geom_line(data = my.rolling.avgs.df, aes(x = Date.Time, y = aou_corrected), color = "black", linewidth = 1) +
  # geom_line(data = my.rolling.avgs.df, aes(x = Date.Time, y = aop), color = col1) +
  # geom_line(data = my.rolling.avgs.df, aes(x = Date.Time, y = o2_bio), color = col3) +
  geom_hline(aes(yintercept = 0), alpha = 0.2) +
  labs(x = "Date", y = expression("Corrected AOP  ["*mu*"M]")) +
  # ylim(c(-100,100)) +
  theme_bw() +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  geom_jitter(aes(x = parse_date_time(sccoos.dates, orders = "Ymd"), y = 140), height = 10, alpha = 0.7) +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12)) 



# ---- cartoon model ----

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

saveRDS(full.aop.df, "full_aop_df.rds")
saveRDS(full.predictors, "full_predictors_df.rds")
saveRDS(mims.hourly, "mims_hourly.rds")


#### to do items ####

## What drive the difference between aop and aop_corrected?  This can be explored in a number of ways, ideally through an ordination (e.g. PCA)
## of the predictors following proper normalization.  Then, look to see if aop-aop_corrected correlates to any of the principal components.