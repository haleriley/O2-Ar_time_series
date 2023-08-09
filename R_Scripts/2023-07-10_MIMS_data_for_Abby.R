
getwd()
setwd('C://Users/haler/Documents/PhD-Bowman/MIMS_O2-Ar_time_series/R_Data/') # set to whatever wd you want
set.seed(1234) # makes random results reproduceable   


# ---- library ----

# install.packages("library.name")

library(tidyverse)
library(lubridate)
library(zoo)
# library(aod)
# library(ranger)
# library(gbm)


# col1.aou <- "#648fff"
# col2.o2bio <- "#dc267f"
# col3.aoucor <- "#fe6100"
# col4.aoupred <- "#785ef0"
# col5.other <- "#ffb000"


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
mims2 <- read.csv('o2bio_vol2.csv')
mims <- rbind(mims, mims2)
mims$date_time <- strptime(mims$date_time, format = '%Y-%m-%d %H')
mims1 <- mims[mims$N2.Ar < 40 & mims$N2.Ar > 30 & mims$date_time < strptime('2021-3-26', format = '%Y-%m-%d'),]
mims2 <- mims[mims$N2.Ar < 20 & mims$N2.Ar > 9 & mims$date_time >= strptime('2021-3-26', format = '%Y-%m-%d'),]
mims <- rbind(mims1, mims2)
mims.date <- mims$date_time
# mims$date_time <- NULL

mims.hourly <- mims %>% group_by(date_time) %>% 
  summarize(temperature = mean(temperature), O2_sat = mean(O2_sat), Ar_sat = mean(Ar_sat), O2.Ar_sat = mean(O2.Ar_sat), 
            O2 = mean(O2), O2.Ar = mean(O2.Ar), N2.Ar = mean(N2.Ar), O2_CF = mean(O2_CF), o2_bio = mean(o2_bio))

plot(mims.date, mims$o2_bio)
# points(mims.hourly.date, mims.hourly$o2_bio,
#        type = 'l',
#        col = 'red')

plot(mims.date, mims$N2.Ar)



### Get miniDOT data ###

miniDot <- read.csv("SIOPier_miniDOT_20180809_20230131_JBowman.csv", header = F, na.strings = "NA", skip = 2)
miniDot <- miniDot[,-1]
colnames(miniDot) <- c("UTC_Date_._Time", "Pacific.Standard.Time", "Temperature", "Dissolved.Oxygen", "Dissolved.Oxygen.Saturation", "Sensor")
# miniDot.colo <- read.csv('../SIOPier_miniDOT_20210430_20210709_BowmanLabSensor.csv', header = T, na.string = "NA")

miniDot$Pacific.Standard.Time <- strptime(miniDot$Pacific.Standard.Time, format = '%Y-%m-%d %H')

## MiniDot reports in mg/L, need uMol

miniDot$Dissolved.Oxygen <- (miniDot$Dissolved.Oxygen / (15.999 * 2 * 1000)) * 10 ** 6

miniDot.col.select <- c("Dissolved.Oxygen", "Temperature")

miniDot.hourly <- miniDot %>% group_by(Pacific.Standard.Time) %>% 
  summarize(Temperature = mean(Temperature), Dissolved.Oxygen = mean(Dissolved.Oxygen), Dissolved.Oxygen.Saturation = mean(Dissolved.Oxygen.Saturation))


### sccoos data ###

## this command should work but doesn't, download with wget instead
sccoos <- read.csv('https://erddap.sccoos.org/erddap/tabledap/autoss.csv?station%2Ctime%2Ctemperature%2Cconductivity%2Cpressure%2Cchlorophyll%2Csalinity&station=%22scripps_pier%22')
# sccoos <- read.csv(file = "autoss_cd8b_3ba9_f9f5.csv")
sccoos <- sccoos[-1,]
sccoos.columns <- c('station', 'time', 'temperature', 'conductivity','pressure','chlorophyll','salinity')
# sccoos <- read.csv('sccoos.csv', header = T, skip = 2, col.names = sccoos.columns)
sccoos$station <- NULL
sccoos$time <- strptime(sccoos$time, format = '%Y-%m-%dT%H', tz = 'GMT')
sccoos$time <- as.POSIXlt(sccoos$time, tz = 'PST')
# sccoos.col.select <- c('temperature', 'pressure','chlorophyll','salinity')
sccoos$temperature <- as.numeric(sccoos$temperature)
sccoos$pressure <- as.numeric(sccoos$pressure)
sccoos$chlorophyll <- as.numeric(sccoos$chlorophyll)
sccoos$conductivity <- as.numeric(sccoos$conductivity)
sccoos.hourly <- sccoos %>% group_by(time) %>% 
  summarize(temperature.sccoos = mean(temperature), conductivity = mean(conductivity), pressure = mean(pressure), chlorophyll = mean(chlorophyll), salinity = mean(salinity))

# sccoos.daily <- sccoos
# sccoos.daily$sample.date <- parse_date_time(paste(year(sccoos.daily$time), month(sccoos.daily$time), day(sccoos.daily$time), sep = "-"), orders = "Ymd")
# sccoos.daily <- sccoos.daily %>% group_by(sample.date) %>%
#   summarize(temperature.sccoos = mean(temperature), conductivity = mean(conductivity), pressure = mean(pressure), chlorophyll = mean(chlorophyll), salinity = mean(salinity))
# saveRDS(sccoos.daily, file = "2023-04-12_sccoos_env_data_daily_mean.rds")



### NOAA - ljac ###

noaa.columns <- c('YY','MM','DD','hh','mm','WDIR','WSPD','GST','WVHT','DPD','APD','MWD','PRES','ATMP','WTMP','DEWP','VIS','TIDE')
noaa.na <- c('99.00', '999', '999.0', '9999.0', '99.0')

ljac.col.select <- c('WSPD','GST','PRES','ATMP','WTMP')
# ljac <- read.table('NOAA_ljac/combined_ljac.txt', col.names = noaa.columns, na.strings = noaa.na)
file.names <- list.files(path = "NOAA_LJAC1_RJH/")
for (i in seq_along(file.names)) {
  assign(substr(file.names[i], start = 1, stop = nchar(file.names[1])-4), read.table(paste("NOAA_LJAC1_RJH/", file.names[i], sep = ""), col.names = noaa.columns, na.strings = noaa.na))
}
ljac <- rbind(ljac1h2018, ljac1h2019, ljac1h2020, ljac1h2021, ljac1h2022, ljac1j2023, ljac1f2023)
ljac$date.time <- paste0(ljac$YY, '-', ljac$MM, '-', ljac$DD, ' ', ljac$hh)
ljac$date.time <- strptime(ljac$date.time, format = '%Y-%m-%d %H', tz = 'GMT')
ljac$date.time <- as.POSIXlt(ljac$date.time, tz = 'PST')
ljac.hourly <- ljac %>% group_by(date.time) %>% 
  summarize(WSPD = mean(WSPD), GST = mean(GST), PRES = mean(PRES), ATMP = mean(ATMP), WTMP = mean(WTMP))

### NOAA - ljpc1 ###

## https://www.ndbc.noaa.gov/station_history.php?station=ljpc1

ljpc1.col.select <- c('WVHT','DPD')
file.names <- list.files(path = "NOAA_LJPC1_RJH/")
for (i in seq_along(file.names)) {
  assign(substr(file.names[i], start = 1, stop = nchar(file.names[1])-4), read.table(paste("NOAA_LJPC1_RJH/", file.names[i], sep = ""), col.names = noaa.columns, na.strings = noaa.na))
}
ljpc1 <- rbind(ljpc1h2018, ljpc1h2019, ljpc1h2020, ljpc1h2021, ljpc1h2022, ljpc1j2023, ljpc1f2023)
ljpc1$date.time <- paste0(ljpc1$YY, '-', ljpc1$MM, '-', ljpc1$DD, ' ', ljpc1$hh)
ljpc1$date.time <- strptime(ljpc1$date.time, format = '%Y-%m-%d %H', tz = 'GMT')
ljpc1$date.time <- as.POSIXlt(ljpc1$date.time, tz = 'PST')
ljpc1.hourly <- ljpc1 %>% group_by(date.time) %>% 
  summarize(WVHT = mean(WVHT), DPD = mean(DPD))


#### combine datasets ####

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

combined.df$Date.Time <- parse_date_time(combined.df$Date.Time, orders = "Ymd HMS")


## The aou calculation is done here, because you either need the O2_sat value
## that comes with the MIMS data or salinity which comes from SCCOOS

combined.df$O2_sat <- O2sat(t = combined.df$Temperature)
combined.df$aou <- -1 * (combined.df$O2_sat - combined.df$Dissolved.Oxygen)
combined.df$delta <- combined.df$o2_bio - combined.df$aou
# really only important for the model part of the project

### some quality checks, could probably do this earlier ###

hist(combined.df$ATMP)
hist(combined.df$PRES)
hist(combined.df$DPD)
hist(combined.df$WVHT)
hist(combined.df$GST)
hist(combined.df$WSPD)
hist(combined.df$pressure)
hist(combined.df$salinity, breaks = 100)
# hist(combined.df$temperature.sccoos.zoo)
hist(combined.df$Dissolved.Oxygen, breaks = 1000)

plot(combined.df$Date.Time, combined.df$Dissolved.Oxygen, type = 'l')

plot(combined.df$o2_bio ~ combined.df$aou)
plot(combined.df$Date.Time, combined.df$aou)

plot(combined.df$Date.Time, combined.df$aou)

points(combined.df$Date.Time, combined.df$o2_bio,
       col = 'red')

## salinity has many problematic values, exclude bad time points
plot(combined.df$salinity~as.Date(combined.df$Date.Time))

combined.df <- combined.df[which(combined.df$salinity > 33.3 & combined.df$salinity < 33.9),]

plot(combined.df$salinity~as.Date(combined.df$Date.Time))

saveRDS(combined.df, file = "2023-04-13_combined_env_data_hourly.rds")


## aggregate by day

library(lubridate)
library(tidyverse)
combined.df$Date <- parse_date_time(paste(day(combined.df$Date.Time), month(combined.df$Date.Time), year(combined.df$Date.Time), sep = "-"), orders = "dmY")
my.colnames <- colnames(combined.df)[c(2:27)]
# my.cols <- c(1:27, 29:30)
combined.df.daily <- combined.df %>% group_by(Date) %>% summarize_at(vars(my.colnames), mean)
saveRDS(combined.df.daily, file = "2023-05-11_combined_env_data_daily.rds")


## nice plots of available data

plot(combined.df$Date.Time, combined.df$aou,
     type = 'l',
     ylab = expression("[O"[2]*"]"[bio]*"|AOU  ["*mu*"M]"),
     xlab = 'Date')

points(combined.df$Date.Time, combined.df$o2_bio, type = 'l', col = 'blue')

legend('topleft',
       legend = c('AOU', expression("[O"[2]*"]"[bio])),
       lty = 1,
       col = c('black', 'blue'))


ggplot(data = combined.df) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  geom_line(aes(x = Date.Time, y = aou), color = "red", lwd = 1, alpha = 0.7) +
  geom_line(aes(x = Date.Time, y = o2_bio), color = "blue", lwd = 1, alpha = 0.7) +
  labs(x = "Date", y = expression("Oxygen Anomaly  [ "*mu*"M]")) +
  theme_bw() +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  theme(panel.grid = element_blank()) +
  ylim(c(-250,250)) +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12))


