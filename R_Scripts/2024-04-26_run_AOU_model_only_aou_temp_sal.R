# BOU estimation GBM, v11
# 2024-08-26
# RJH

## MIMS v2.1 file download changed, so making a new script 


## A more correct approach would be to derive a "true" O2 value for the MIMS
## based on existing calibrations, then base the model only on timepoints where
## the mims and miniDOT O2 values agree
getwd()
setwd('C://Users/haler/Documents/PhD-Bowman/O2-Ar_time_series/R_Data/')
set.seed(1234)

predictors <- c('aou',
                #'delta',
                "temperature.sccoos",
                'salinity',
                # "pressure",
                # "WSPD",
                # "GST",
                # "WVHT",
                # "DPD",
                ##"ATMP",
                ##"PRES",
                "BOU")

accessory <- c('Date.Time', 'Dissolved.Oxygen', 'O2')



# ---- library ----

library(tidyverse)
library(lubridate)
library(zoo)
library(aod)
library(ranger)
library(gbm)
library(rLakeAnalyzer)
library(readxl)
library(vegan)
# library(plotly)
library(patchwork)
library(RColorBrewer)


col1.aou <- "#648fff"
col2.o2bio <- "#dc267f"
col3.BOU_est <- "#fe6100" 
col4.BOU_pred <- "#785ef0" 
col5.other <- "#ffb000"



# ---- functions ----

## Define function to calculate Ar at saturation based on Hamme and Emerson, 2004
Arsat <- function(salinity,temperature){
  
  TS = log((298.15-temperature) / (273.15+temperature))
  
  A0 = 2.7915
  A1 = 3.17609
  A2 = 4.13116
  A3 = 4.90379
  B0 = -6.96233 * 10 ** -3
  B1 = -7.66670 * 10 ** -3
  B2 = -1.16888 * 10 ** -2
  
  Ar = exp(A0 + A1*TS + A2*TS^2 + A3*TS^3 + salinity*(B0 + B1*TS + B2*TS^2))
  
  ## final units are umol kg-1
  
  return(Ar)
  
}

## Define function to calculate O2 at saturation based on Garcia and Gordon, 1992
O2sat <- function(salinity, temperature){
  
  TS = log((298.15-temperature) / (273.15+temperature))
  
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
  
  O2 = exp(A0 + A1*TS + A2*TS^2 + A3*TS^3 + A4*TS^4 + A5*TS^5 + salinity*(B0 + B1*TS + B2*TS^2 + B3*TS^3 + C0*salinity^2))
  
  ## final units are umol kg-1
  
  return(O2)
  
}

set.plot.groups <- function(df){
  
  df$time.change <- difftime(c(df$Date.Time[-1], NA), df$Date.Time)
  index <- c(1,1+which(df$time.change > 24))
  g <- 1
  for(i in index){
    
    df$plot.group[i:nrow(df)] <- g
    g <- g+1
    
  }
  
  return(df)
}


# ---- ignore for this one ----
#' # ---- testing median salinity on [O2]sat ----
#' 
#' ## test if median salinity is okay to use for calculating [O2]sat
#' combined.df <- readRDS("2024-01-19_combined_env_data_hourly.rds")
#' test.temps <- range(na.omit(combined.df$temperature))
#' test.lowS <- O2sat(33.0, test.temps)
#' test.highS <- O2sat(33.9, test.temps)
#' 
#' 
#' 
#' # ---- load data ----
#' 
#' ### Get miniDOT data ###
#' 
#' miniDot <- read.csv("Current_Model_Inputs/miniDOT/SIOPier_miniDOT_20180809_20230503_JBowman (2).csv", header = F, na.strings = "NA", skip = 2)
#' miniDot2 <- read_excel(path = "Current_Model_Inputs/miniDOT/SIOPier_miniDOT_20230503_20231005_BowmanLab.xlsx")
#' miniDot2 <- miniDot2[-1,-c(4,8)]
#' miniDot3 <- read.csv("Current_Model_Inputs/miniDOT/SIOPier_miniDOT_20231005_20240524_Timeseries_RHale.csv", header = F, na.strings = "NA", skip = 1)
#' 
#' miniDot <- miniDot[,-1]
#' miniDot2 <- miniDot2[,-1]
#' miniDot3 <- miniDot3[,-1]
#' 
#' colnames(miniDot) <- c("UTC_Date_._Time", "Pacific.Standard.Time", "Temperature", "Dissolved.Oxygen", "Dissolved.Oxygen.Saturation", "Sensor")
#' colnames(miniDot2) <- c("UTC_Date_._Time", "Pacific.Standard.Time", "Temperature", "Dissolved.Oxygen", "Dissolved.Oxygen.Saturation", "Sensor")
#' colnames(miniDot3) <- c("UTC_Date_._Time", "Pacific.Standard.Time", "Temperature", "Dissolved.Oxygen", "Dissolved.Oxygen.Saturation", "Sensor")
#' 
#' 
#' miniDot <- rbind(miniDot, miniDot2, miniDot3)
#' miniDot$Temperature <- as.numeric(miniDot$Temperature)
#' miniDot$Dissolved.Oxygen <- as.numeric(miniDot$Dissolved.Oxygen)
#' miniDot$Dissolved.Oxygen.Saturation <- as.numeric(miniDot$Dissolved.Oxygen.Saturation)
#' 
#' 
#' miniDot$Pacific.Standard.Time <- strptime(miniDot$Pacific.Standard.Time, format = '%Y-%m-%d %H')
#' ## I don't think we need this based off the info Samantha gave me (deployment ends 9pm on 1/31/23) 
#' 
#' ## MiniDot reports in mg/L, need uMol
#' miniDot$Dissolved.Oxygen <- (miniDot$Dissolved.Oxygen / (15.999 * 2 * 1000)) * 10 ** 6
#' 
#' miniDot.col.select <- c("Dissolved.Oxygen", "Temperature")
#' 
#' miniDot <- na.omit(miniDot)
#' 
#' ## aggregate hourly ##
#' miniDot.hourly <- miniDot %>% group_by(Pacific.Standard.Time) %>% 
#'   summarize(Temperature = mean(Temperature), Dissolved.Oxygen = mean(Dissolved.Oxygen), Dissolved.Oxygen.Saturation = mean(Dissolved.Oxygen.Saturation))
#' 
#' 
#' 
#' ### Get sccoos data ###
#' # sccoos1<- read.csv('https://erddap.sccoos.org/erddap/tabledap/autoss.csv?station%2Ctime%2Ctemperature%2Cconductivity%2Cpressure%2Cchlorophyll%2Csalinity&station=%22scripps_pier%22&time%3E=%202005-06-16T19%3A27%3A53Z&time%3C=2010-12-10T00%3A00%3A00Z') # file will not download if it is too big
#' # sccoos2<- read.csv('https://erddap.sccoos.org/erddap/tabledap/autoss.csv?station%2Ctime%2Ctemperature%2Cconductivity%2Cpressure%2Cchlorophyll%2Csalinity&station=%22scripps_pier%22&time%3E=2010-12-10T00%3A00%3A00Z&time%3C=2017-12-10T00%3A00%3A00Z') # file will not download if it is too big
#' # sccoos3<- read.csv('https://erddap.sccoos.org/erddap/tabledap/autoss.csv?station%2Ctime%2Ctemperature%2Cconductivity%2Cpressure%2Cchlorophyll%2Csalinity&station=%22scripps_pier%22&time%3E=2017-12-10T00%3A00%3A00Z&time%3C=2024-06-27T05%3A08%3A00Z') # file will not download if it is too big
#' # sccoos <- rbind(sccoos1, sccoos2, sccoos3)
#' # sccoos <- read.csv('https://erddap.sccoos.org/erddap/tabledap/autoss.csv?station%2Ctime%2Ctemperature%2Cconductivity%2Cpressure%2Cchlorophyll%2Csalinity&station=%22scripps_pier%22')
#' 
#' options(timeout = 300)
#' sccoos <- read.csv("https://erddap.sensors.axds.co/erddap/tabledap/scripps-pier-automated-shore-sta-1.csv?time%2Cmass_concentration_of_chlorophyll_in_sea_water%2Csea_water_electrical_conductivity%2Csea_water_practical_salinity%2Csea_water_pressure%2Csea_water_temperature%2Cstation")
#' 
#' sccoos <- sccoos[-1,]
#' colnames(sccoos) <- c('time', 'chlorophyll', 'conductivity', 'salinity','pressure',  'temperature', 'station')
#' sccoos$station <- NULL
#' sccoos$time <- strptime(sccoos$time, format = '%Y-%m-%dT%H', tz = 'GMT')
#' sccoos$time <- as.POSIXlt(sccoos$time, tz = 'PST')
#' sccoos$temperature <- as.numeric(sccoos$temperature)
#' sccoos$pressure <- as.numeric(sccoos$pressure)
#' sccoos$chlorophyll <- as.numeric(sccoos$chlorophyll)
#' sccoos$conductivity <- as.numeric(sccoos$conductivity)
#' 
#' sccoos <- sccoos[-which(is.nan(sccoos$temperature)),]
#' 
#' ## aggregate hourly ##
#' sccoos.hourly <- sccoos %>% group_by(time) %>% 
#'   summarize(temperature.sccoos = mean(temperature), conductivity = mean(conductivity), pressure = mean(pressure), chlorophyll = mean(chlorophyll), salinity = mean(salinity))
#' 
#' 
#' ## can aggregate daily if desired ##
#' 
#' # sccoos.daily <- sccoos
#' # sccoos.daily$sample.date <- parse_date_time(paste(year(sccoos.daily$time), month(sccoos.daily$time), day(sccoos.daily$time), sep = "-"), orders = "Ymd")
#' # sccoos.daily <- sccoos.daily %>% group_by(sample.date) %>%
#' #   summarize(temperature.sccoos = mean(temperature), conductivity = mean(conductivity), pressure = mean(pressure), chlorophyll = mean(chlorophyll), salinity = mean(salinity))
#' # saveRDS(sccoos.daily, file = "2024-02-14_sccoos_env_data_daily_mean.rds")
#' 
#' 
#' ## Get O2bio calculated in read_lvm.py ##
#' 
#' mims <- read.csv('Current_Model_Inputs/MIMS_o2bio/o2bio.csv') # this file is produced by read_lvm.py, refer to that script for details
#' 
#' mims$salinity.MIMS <- NA
#' mims$O2.MIMS <- NA
#' mims$O2_sat.MIMS <- NA
#' mims$fluorescence.MIMS <- NA
#' colnames(mims)[which(colnames(mims) == "temperature")] <- "Inlet.Temperature"
#' colnames(mims)[which(colnames(mims) == "o2_bio")] <- "BOU"
#' 
#' # temp <- tempfile()
#' # download.file("https://polarmicrobes.org/ecoobs/MIMS_data_vol2.csv.gz", temp) # loads in current MIMS vol2 data
#' # my.file <- read.csv(gzfile(temp), as.is = TRUE)
#' # unlink(temp)
#' # my.file <- read.csv("2024-02-16_MIMS_data_vol2.csv") # THIS NEEDS TO BE FIXED!!
#' 
#' temp <- tempfile()
#' download.file("https://polarmicrobes.org/ecoobs/o2bio_vol2.csv", temp) # loads in current MIMS vol2 data
#' my.file <- read.csv(gzfile(temp), as.is = TRUE)
#' unlink(temp)
#' 
#' mims2 <- my.file
#' # mims2 <- my.file[,c("date_time", "O2", "Temperature..ITS.90.deg.C.")]
#' mims2$date_time <- parse_date_time(paste(year(mims2$date_time), "-", month(mims2$date_time), "-", day(mims2$date_time), " ", hour(mims2$date_time), sep = ""), orders = "Ymd H")
#' # mims2 <- my.file[,c("date_time", "N2", "O2", "Ar", "Inlet.Temperature")]
#' # mims2$date_time <- parse_date_time(paste(year(mims2$time), "-", month(mims2$time), "-", day(mims2$time), " ", hour(mims2$time), sep = ""), orders = "Ymd H")
#' # mims2 <- mims2[,-1]
#' mims2 <- mims2 %>% group_by(date_time) %>% summarize_all(mean)
#' colnames(mims2)[which(colnames(mims2) == "o2_bio")] <- "BOU"
#' mims2$salinity.MIMS <- NA
#' mims2$O2.MIMS <- NA
#' mims2$O2_sat.MIMS <- NA
#' mims2$fluorescence.MIMS <- NA
#' colnames(mims2)[which(colnames(mims2) == "temperature")] <- "Inlet.Temperature"
#' 
#' 
#' 
#' temp <- tempfile()
#' download.file("https://polarmicrobes.org/ecoobs/o2bio_vol2.1.csv", temp) # loads in current MIMS vol2 data
#' my.file <- read.csv(gzfile(temp), as.is = TRUE)
#' unlink(temp)
#' 
#' mims2.1 <- my.file
#' # mims2 <- my.file[,c("date_time", "O2", "Temperature..ITS.90.deg.C.")]
#' mims2.1$date_time <- parse_date_time(paste(year(mims2.1$date_time), "-", month(mims2.1$date_time), "-", day(mims2.1$date_time), " ", hour(mims2.1$date_time), sep = ""), orders = "Ymd H")
#' # mims2 <- my.file[,c("date_time", "N2", "O2", "Ar", "Inlet.Temperature")]
#' # mims2$date_time <- parse_date_time(paste(year(mims2$time), "-", month(mims2$time), "-", day(mims2$time), " ", hour(mims2$time), sep = ""), orders = "Ymd H")
#' # mims2 <- mims2[,-1]
#' mims2.1 <- mims2.1 %>% group_by(date_time) %>% summarize_all(mean)
#' 
#' colnames(mims2.1)[2:6] <- c("Inlet.Temperature", "salinity.MIMS", "O2.MIMS", "O2_sat.MIMS", "fluorescence.MIMS")
#' colnames(mims2.1)[which(colnames(mims2.1) == "o2_bio")] <- "BOU"
#' 
#' 
#' # temp.sccoos <- sccoos[,c("time", "temperature", "salinity")]
#' # colnames(temp.sccoos)[1] <- "date_time"
#' # temp.sccoos$date_time <- parse_date_time(as.character(temp.sccoos$date_time), orders = "Ymd HMS")
#' # temp.sccoos <- temp.sccoos %>% group_by(date_time) %>% summarize_all(mean)
#' # mims2.merged <- merge(mims2, temp.sccoos, by = "date_time", all.x = T, all.y = T)
#' 
#' ## choosing to use SCCOOS data because they are collected together
#' # mims2.merged$Ar_sat <- Arsat(salinity = mims2.merged$salinity, temperature = mims2.merged$temperature)
#' # mims2.merged$O2_sat <- O2sat(salinity = mims2.merged$salinity, temperature = mims2.merged$temperature)  
#' # mims2.merged$O2.Ar_sat <- mims2.merged$O2_sat/mims2.merged$Ar_sat
#' # mims2.merged$O2.Ar <- mims2.merged$O2/mims2.merged$Ar
#' # mims2.merged$N2.Ar <- mims2.merged$N2/mims2.merged$Ar
#' # mims2.merged$O2_CF <- 1.54
#' # mims2.merged$O2_CF[which(mims2.merged$date_time >= parse_date_time('2022-11-22 12:00:00', orders = "Ymd HMS"))] <- 1.76
#' # mims2.merged$O2_CF[which(mims2.merged$date_time >= parse_date_time('2023-01-23 12:00:00', orders = "Ymd HMS"))] <- 2.0
#' # mims2.merged$O2_CF[which(mims2.merged$date_time >= parse_date_time('2023-03-29 12:00:00', orders = "Ymd HMS"))] <- 1.83
#' # mims2.merged$O2_CF[which(mims2.merged$date_time >= parse_date_time('2023-05-12 12:00:00', orders = "Ymd HMS"))] <- 1.68
#' # mims2.merged$O2_CF[which(mims2.merged$date_time >= parse_date_time('2023-10-09 12:00:00', orders = "Ymd HMS"))] <- 1.54
#' # 
#' # mims2.merged$BOU <- ((((mims2.merged$O2.Ar * mims2.merged$O2_CF) / mims2.merged$O2.Ar_sat) - 1) * mims2.merged$O2_sat * 1) # assuming [Ar]/[Ar]sat = 1, from Cassar et al., 2011
#' # ## postitive BOU here = net autotrophy/supersaturation, same as MIMS1/read_lvm.py script
#' # mims2 <- mims2.merged[,-c(6:7)]
#' 
#' # ## for Abby
#' # mims2$date <- parse_date_time(substr(as.character(mims2$date_time), start = 1, stop = 10), orders = "Ymd")
#' # temp <- mims2 %>% group_by(date) %>% summarize_all(mean)
#' # write.csv(temp, "2024-06-26_O2_Ar_daily_for_Abby.csv")
#' 
#' 
#' mims <- mims[,colnames(mims2.1)]
#' mims2 <- mims2[,colnames(mims2.1)]
#' 
#' 
#' mims$date_time <- parse_date_time(paste(year(mims$date_time), "-", month(mims$date_time), "-", day(mims$date_time), " ", hour(mims$date_time), sep = ""), orders = "Ymd H")
#' mims <- mims %>% group_by(date_time) %>% summarize_all(mean)
#' 
#' # sat.df <- mims2.merged[,c("date_time", "Ar_sat", "O2_sat", "O2.Ar_sat")]
#' 
#' mims <- rbind(mims, mims2, mims2.1) # combine datasets
#' # mims$date_time <- strptime(mims$date_time, format = '%Y-%m-%d %H')
#' mimsQC1 <- mims[mims$N2.Ar < 40 & mims$N2.Ar > 30 & mims$date_time < parse_date_time('2021-3-26 0-0-0', orders = "Ymd HMS"),] # manual QC
#' mimsQC2 <- mims[mims$N2.Ar < 20 & mims$N2.Ar > 9 & mims$date_time >= parse_date_time('2021-3-26 0-0-0', orders = "Ymd HMS"),]
#' mims <- rbind(mimsQC1, mimsQC2)
#' mims.date <- mims$date_time
#' # mims <- mims[,-which(colnames(mims) %in% c("Ar_sat", "O2_sat", "O2.Ar_sat"))]
#' 
#' 
#' 
#' 
#' mims.hourly <- mims ## just renaming 
#' 
#' ## quick data visualization
#' plot(mims.date, mims$BOU)
#' plot(mims.date, mims$N2.Ar)
#' 
#' 
#' 
#' ### Get NOAA LJAC1 buoy data ###
#' ## Data downloaded from: https://www.ndbc.noaa.gov/station_history.php?station=ljac1
#' 
#' noaa.columns <- c('YY','MM','DD','hh','mm','WDIR','WSPD','GST','WVHT','DPD','APD','MWD','PRES','ATMP','WTMP','DEWP','VIS','TIDE')
#' noaa.na <- c('99.00', '999', '999.0', '9999.0', '99.0')
#' 
#' ljac.col.select <- c('WSPD','GST','PRES','ATMP','WTMP')
#' 
#' # files all saved in same folder
#' file.names <- list.files(path = "Current_Model_Inputs/NOAA_LJAC1_RJH/")
#' for (i in seq_along(file.names)) {
#'   assign(substr(file.names[i], start = 1, stop = nchar(file.names[i])-4), read.table(paste("Current_Model_Inputs/NOAA_LJAC1_RJH/", file.names[i], sep = ""), col.names = noaa.columns, na.strings = noaa.na))
#' }
#' 
#' ljac <- rbind(ljac1h2018, ljac1h2019, ljac1h2020, ljac1h2021, ljac1h2022, ljac1h2023, ljac1jan2024, ljac1feb2024, ljac1mar2024, ljac1apr2024, ljac1may2024, ljac1jun2024)
#' ljac$date.time <- paste0(ljac$YY, '-', ljac$MM, '-', ljac$DD, ' ', ljac$hh)
#' ljac$date.time <- strptime(ljac$date.time, format = '%Y-%m-%d %H', tz = 'GMT')
#' ljac$date.time <- as.POSIXlt(ljac$date.time, tz = 'PST')
#' 
#' ## aggregate hourly
#' ljac.hourly <- ljac %>% group_by(date.time) %>% 
#'   summarize(WSPD = mean(WSPD), GST = mean(GST), PRES = mean(PRES), ATMP = mean(ATMP), WTMP = mean(WTMP))
#' 
#' 
#' ### Get NOAA LJPC1 buoy data ###
#' ## Data downloaded from: https://www.ndbc.noaa.gov/station_history.php?station=ljpc1
#' 
#' noaa.columns <- c('YY','MM','DD','hh','mm','WDIR','WSPD','GST','WVHT','DPD','APD','MWD','PRES','ATMP','WTMP','DEWP','VIS','TIDE')
#' noaa.na <- c('99.00', '999', '999.0', '9999.0', '99.0')
#' 
#' ljpc1.col.select <- c('WVHT','DPD')
#' 
#' # files all saved in same folder
#' file.names <- list.files(path = "Current_Model_Inputs/NOAA_LJPC1_RJH/")
#' for (i in seq_along(file.names)) {
#'   assign(substr(file.names[i], start = 1, stop = nchar(file.names[i])-4), read.table(paste("Current_Model_Inputs/NOAA_LJPC1_RJH/", file.names[i], sep = ""), col.names = noaa.columns, na.strings = noaa.na))
#' }
#' 
#' ljpc1 <- rbind(ljpc1h2018, ljpc1h2019, ljpc1h2020, ljpc1h2021, ljpc1h2022, ljpc1h2023, ljpc1jan2024, ljpc1feb2024, ljpc1mar2024, ljpc1apr2024, ljpc1may2024, ljpc1jun2024)
#' ljpc1$date.time <- paste0(ljpc1$YY, '-', ljpc1$MM, '-', ljpc1$DD, ' ', ljpc1$hh)
#' ljpc1$date.time <- strptime(ljpc1$date.time, format = '%Y-%m-%d %H', tz = 'GMT')
#' ljpc1$date.time <- as.POSIXlt(ljpc1$date.time, tz = 'PST')
#' 
#' ## aggregate hourly
#' ljpc1.hourly <- ljpc1 %>% group_by(date.time) %>% 
#'   summarize(WVHT = mean(WVHT), DPD = mean(DPD))
#' 
#' 
#' 
#' # ---- combine datasets ----
#' 
#' colnames(mims.hourly)[1] <- "Date.Time"
#' colnames(miniDot.hourly)[1] <- "Date.Time"
#' colnames(sccoos.hourly)[1] <- "Date.Time"
#' colnames(ljac.hourly)[1] <- "Date.Time"
#' colnames(ljpc1.hourly)[1] <- "Date.Time"
#' # colnames(sat.df)[1] <- "Date.Time"
#' 
#' mims.hourly <- data.frame(mims.hourly)
#' miniDot.hourly <- data.frame(miniDot.hourly)
#' sccoos.hourly <- data.frame(sccoos.hourly)
#' ljac.hourly <- data.frame(ljac.hourly)
#' ljpc1.hourly <- data.frame(ljpc1.hourly)
#' # sat.df <- data.frame(sat.df)
#' 
#' mims.hourly$Date.Time <- as.character(mims.hourly$Date.Time)
#' miniDot.hourly$Date.Time <- as.character(miniDot.hourly$Date.Time)
#' sccoos.hourly$Date.Time <- as.character(sccoos.hourly$Date.Time)
#' ljac.hourly$Date.Time <- as.character(ljac.hourly$Date.Time)
#' ljpc1.hourly$Date.Time <- as.character(ljpc1.hourly$Date.Time)
#' # sat.df$Date.Time <- as.character(sat.df$Date.Time)
#' 
#' combined.df <- merge(mims.hourly, miniDot.hourly, by = "Date.Time", all = T, sort = T)
#' combined.df <- merge(combined.df, sccoos.hourly, by = "Date.Time", all = T, sort = T)
#' combined.df <- merge(combined.df, ljac.hourly, by = "Date.Time", all = T, sort = T)
#' combined.df <- merge(combined.df, ljpc1.hourly, by = "Date.Time", all = T, sort = T)
#' # combined.df <- merge(combined.df, sat.df, by = "Date.Time", all = T, sort = T)
#' 
#' ## turn date into PociXct format using parse_date_time() function
#' combined.df$Date.Time <- parse_date_time(combined.df$Date.Time, orders = "Ymd HMS")
#' 
#' combined.df <- combined.df[which(is.na(combined.df$Dissolved.Oxygen) == F),]
#' 
#' combined.df$O2_sat <- O2sat(salinity = combined.df$salinity, temperature = combined.df$temperature.sccoos)
#' 
#' 
#' # ---- calculate aou and delta ----
#' 
#' ## The aou calculation is done here, because you either need the O2_sat value
#' ## that comes with the MIMS data or salinity which comes from SCCOOS
#' 
#' combined.df$aou <- -(combined.df$O2_sat - combined.df$Dissolved.Oxygen)
#' combined.df$delta <- combined.df$BOU - combined.df$aou # hmm but BOU is also calculated with O2_sat in the python script
#' 
#' 
#' # ---- QC data ----
#' 
#' hist(combined.df$ATMP)
#' hist(combined.df$PRES)
#' hist(combined.df$DPD)
#' hist(combined.df$WVHT)
#' hist(combined.df$GST)
#' hist(combined.df$WSPD)
#' hist(combined.df$pressure)
#' hist(combined.df$salinity, breaks = 100)
#' hist(combined.df$temperature.sccoos)
#' hist(combined.df$Dissolved.Oxygen, breaks = 1000)
#' 
#' plot(combined.df$Date.Time, combined.df$Dissolved.Oxygen, type = 'l')
#' 
#' plot(combined.df$BOU ~ combined.df$aou)
#' plot(combined.df$Date.Time, combined.df$aou)
#' 
#' ## compare MIMS o2bio time series and miniDOT aou time series
#' plot(combined.df$Date.Time, combined.df$aou)
#' points(combined.df$Date.Time, combined.df$BOU,
#'        col = 'red')
#' 
#' 
#' ## salinity has many problematic values, exclude bad time points
#' plot(combined.df$salinity~as.Date(combined.df$Date.Time))
#' combined.df <- combined.df[which(combined.df$salinity > 33 & combined.df$salinity < 33.8),]
#' plot(combined.df$salinity~as.Date(combined.df$Date.Time))
#' 
#' # remove unreasonably high [O2]bio measurements
#' combined.df <- combined.df[which(combined.df$BOU < 400 | is.na(combined.df$BOU)),]
#' 
#' 
#' ## check to see where O2 and O2 agree ##
#' 
#' plot(combined.df$O2[combined.df$Date.Time < strptime('2021-3-26', format = '%Y-%m-%d')],
#'      combined.df$Dissolved.Oxygen[combined.df$Date.Time < strptime('2021-3-26', format = '%Y-%m-%d')],
#'      log = 'x')
#' 
#' plot(combined.df$O2,
#'      combined.df$Dissolved.Oxygen,
#'      log = 'x')
#' 
#' 
#' ## data selection for modeling, constrain by vacuum pressure (MIMS) ##
#' 
#' index <- which(combined.df$O2 < 1.5e-9)
#' combined.df$O2[index] <- NA
#' combined.df$BOU[index] <- NA
#' 
#' predictors <- c('aou',
#'                 #'delta',
#'                 "temperature.sccoos",
#'                 'salinity',
#'                 "pressure",
#'                 "WSPD",
#'                 "GST",
#'                 "WVHT",
#'                 "DPD",
#'                 #"ATMP",
#'                 #"PRES",
#'                 "BOU")
#' 
#' accessory <- c('Date.Time', 'Dissolved.Oxygen', 'O2')
#' 
#' combined.df <- combined.df[,c(predictors, accessory)]


temp.cols <- c(predictors, accessory)
temp.cols <- temp.cols[-which(temp.cols == "BOU" | temp.cols == "O2")]

#' combined.df <- combined.df %>% drop_na(temp.cols) # ignore warning
#' 
#' 
#' # saveRDS(combined.df, file = "2024-08-26_combined_env_data_hourly.rds")
#' 
#' 
#' 
#' # ---- aggregate by day if desired ----
#' 
#' combined.df$Date <- parse_date_time(paste(day(combined.df$Date.Time), month(combined.df$Date.Time), year(combined.df$Date.Time), sep = "-"), orders = "dmY")
#' my.colnames <- colnames(combined.df)[-which(colnames(combined.df) == "Date.Time")]
#' combined.df.daily <- combined.df %>% group_by(Date) %>% summarize_all(mean)
#' # saveRDS(combined.df.daily, file = "2024-08-26_combined_env_data_daily.rds")
#' 

# ---- nice plot of available AOU and BOU data -----


combined.df <- readRDS(file = "2024-08-26_combined_env_data_hourly.rds")

combined.df.overlap <- combined.df[which(is.na(combined.df$BOU) == FALSE),]
combined.df.overlap <- set.plot.groups(combined.df.overlap)

combined.df <- set.plot.groups(combined.df)

# png("../Figures/manuscript_figures/2024-08-26_fig1.png", width = 1000, height = 500)
# # png("../Figures/manuscript_figures/2024-05-17_BOU.png", width = 1000, height = 500)
# 
# # ggplot(data = combined.df.overlap) +
# plot.1a <- ggplot(data = combined.df) +
#   geom_hline(yintercept = 0, alpha = 0.5) +
#   # geom_line(aes(x = Date.Time, y = aou), color = col1.aou, lwd = 1, alpha = 0.7) +
#   geom_line(aes(x = Date.Time, y = aou, color = "AOU", group = plot.group), lwd = 1, alpha = 0.7) +
#   geom_line(aes(x = Date.Time, y = BOU, color = "BOU", group = plot.group), lwd = 1, alpha = 0.7) +
#   labs(x = "Date", y = "Oxygen Anomaly [uM]", color = "") +
#   theme_bw() +
#   scale_x_datetime(date_breaks = "1 year", date_labels = "%Y", 
#                    date_minor_breaks = "3 months",
#                    # limits = c(min(combined.df$Date.Time), max(combined.df$Date.Time))
#   ) +
#   # theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
#   theme(panel.grid = element_blank()) +
#   ylim(c(-250,250)) +
#   scale_color_brewer(palette = "Dark2") +
#   # scale_color_manual(name = "", labels = factor(c("AOU", "BOU"), levels = c("BOU", "AOU")), guide = "legend", values = c("BOU" = col2.o2bio, "AOU" = col1.aou)) +
#   # scale_x_datetime(date_breaks = "3 days", date_labels = "%b %d") +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) 
# plot.1a
# 
# dev.off()
# 
# mean.BOU.full <- mean(combined.df$BOU[which(is.na(combined.df$BOU) == F)])
# sd.BOU.full <- sd(combined.df$BOU[which(is.na(combined.df$BOU) == F)])
# 
# mean.AOU.full <- mean(combined.df$aou)
# sd.AOU.full <- sd(combined.df$aou)


date1 <- parse_date_time('2021-7-1', orders = "Ymd")
date2 <- parse_date_time('2021-7-14', orders = "Ymd")

# plot.2a <- ggplot(data = combined.df.overlap) +
#   # ggplot(data = combined.df) +
#   geom_hline(yintercept = 0, alpha = 0.5) +
#   # geom_line(aes(x = Date.Time, y = aou), color = col1.aou, lwd = 1, alpha = 0.7) +
#   geom_line(aes(x = Date.Time, y = aou, color = "AOU", group = plot.group), lwd = 1, alpha = 0.7) +
#   geom_line(aes(x = Date.Time, y = BOU, color = "BOU", group = plot.group), lwd = 1, alpha = 0.7) +
#   # geom_line(data = combined.df.overlap[which(combined.df.overlap$Date.Time >= date1 & combined.df.overlap$Date.Time <= date2),],
#   #             aes(x = Date.Time, y = BOU), lwd = 1, alpha = 1, color = brewer.pal(8, "Dark2")[6]) +
#   labs(x = "Date", y = "Oxygen Anomaly [uM]", color = "") +
#   theme_bw() +
#   # scale_x_datetime(date_breaks = "1 year", date_labels = "%Y",
#   #                  date_minor_breaks = "3 months",
#   #                  limits = c(min(combined.df$Date.Time), max(combined.df$Date.Time))
#   # ) +
#   # theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
#   theme(panel.grid = element_blank()) +
#   ylim(c(-250,250)) +
#   scale_color_brewer(palette = "Dark2") +
#   # scale_color_manual(name = "", labels = factor(c("AOU", "BOU"), levels = c("BOU", "AOU")), guide = "legend", values = c("BOU" = col2.o2bio, "AOU" = col1.aou)) +
#   # scale_x_datetime(date_breaks = "3 days", date_labels = "%b %d") +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) 
# plot.2a



# ---- establish training and validation data, set up GBM ----

# combined.df <- readRDS(file = "2024-08-26_combined_env_data_hourly.rds")
combined.df.select <- combined.df %>% drop_na(c(temp.cols, "BOU"))
# 
# # saveRDS(combined.df.select, file = "2024-08-26_combined.df.select.rds")
# combined.df.select <- readRDS("2024-08-26_combined.df.select.rds")



## date range for validation, shows dynamics well

combined.df.select$Date.Time <- parse_date_time(combined.df.select$Date.Time, orders = "Ymd HMS")
date1 <- parse_date_time('2021-7-1', orders = "Ymd")
date2 <- parse_date_time('2021-7-14', orders = "Ymd")

combined.df.other <- combined.df.select[which(combined.df.select$Date.Time > parse_date_time('2022-04-15', orders = "Ymd") &
                                                combined.df.select$Date.Time < parse_date_time('2022-04-29', orders = "Ymd")),]
combined.df.other.Date.Time <- combined.df.other$Date.Time


combined.df.test <- combined.df.select[which(combined.df.select$Date.Time > date1 &
                                               combined.df.select$Date.Time < date2),]

combined.df.train <- combined.df.select[which(!combined.df.select$Date.Time %in% combined.df.test$Date.Time),]

combined.df.test.Date.Time <- combined.df.test$Date.Time
combined.df.train.Date.Time <- combined.df.train$Date.Time

combined.df.train <- combined.df.train[,predictors]
combined.df.test <- combined.df.test[,predictors]


m0.rmse.test <- sqrt(mean((combined.df.test$aou - combined.df.test$BOU)^2))
m0.perc.err.test <- median(abs(combined.df.test$aou - combined.df.test$BOU)/abs(combined.df.test$BOU))
m0.R2.test <- summary(lm(combined.df.test$aou ~ combined.df.test$BOU))$r.squared

m0.rmse.full <- sqrt(mean((combined.df.select$aou - combined.df.select$BOU)^2))
m0.perc.err.full <- mean(abs(combined.df.select$aou - combined.df.select$BOU)/abs(combined.df.select$BOU))
m0.R2.full <- summary(lm(combined.df.select$aou ~ combined.df.select$BOU))$r.squared

# # plot all aou and [O2]bio data used in the model
# ggplot() +
#   geom_hline(aes(yintercept = 0), alpha = 0.2) +
#   geom_line(aes(x = combined.df.select$Date.Time[-c(1:2)], y = combined.df.select$BOU[-c(1:2)]), color = col2.o2bio, lwd = 1, alpha = 0.7) +
#   geom_line(aes(x = combined.df.select$Date.Time[-c(1:2)], y = combined.df.select$aou[-c(1:2)]), color = col1.aou, lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df.test.Date.Time, y = combined.df.test$BOU), color = col3, lwd = 1, alpha = 1) +
#   labs(x = "Date", y = expression("[O"[2]*"]"[bio]*" | aou  ["*mu*"M]")) +
#   ylim(c(-250,250)) +
#   theme_bw() +
#   scale_x_datetime(date_breaks = "3 months", date_labels = "%b %Y") +
#   theme(panel.grid = element_blank()) +
#   theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12))


# # ---- PCA of predictors to compare train and test data ----
# 
# my.date.labels <- combined.df.select$Date.Time
# env.pred.matrix <- as.matrix(combined.df.select[,predictors])
# env.pred.matrix <- env.pred.matrix[,-9]
# 
# PCA <- rda(env.pred.matrix, scale = FALSE)
# 
# barplot(as.vector(PCA$CA$eig)/sum(PCA$CA$eig))
# sum((as.vector(PCA$CA$eig)/sum(PCA$CA$eig))[1:2])
# smry <- summary(PCA)
# df1  <- data.frame(smry$sites[,1:2])# PC1 and PC2
# df1$Date.Time <- my.date.labels
# df2  <- data.frame(smry$species[,1:2])     # loadings for PC1 and PC2
# df1$data.set <- "train"
# df1$data.set[which(df1$Date.Time >= date1 & df1$Date.Time <= date2)] <- "test"
# df1$data.set <- factor(df1$data.set, levels = c("train", "test"))
# 
# # df2 <- df2[which((abs(df2$PC1) + abs(df2$PC2)) >= 0.3),]
# 
# plot.3a <- ggplot() + 
#   # geom_point(aes(color = data.set, alpha = data.set)) +
#   geom_point(data = df1[which(df1$data.set == "train"),], aes(x=PC1, y=PC2, color = "Training"), alpha = 0.3) +
#   geom_point(data = df1[which(df1$data.set == "test"),], aes(x=PC1, y=PC2, color = "Validation"), alpha = 0.8) +
#   geom_hline(yintercept=0, linetype="dotted") +
#   geom_vline(xintercept=0, linetype="dotted") +
#   coord_fixed() +
#   # geom_segment(data=df2, aes(x=0, xend=PC1, y=0, yend=PC2),
#   #              color="black", arrow=arrow(length=unit(0.01,"npc"))) +
#   # geom_text(data=df2,
#   #           aes(x=PC1,y=PC2,label=rownames(df2)),
#   #               # hjust= (1*sign(PC1)), vjust= (-0.5*sign(PC1))),
#   #           color="black", size=4) +
#   xlim(c(-max(c(df1$PC1, df1$PC2)), max(c(df1$PC1, df1$PC2)))) + ylim(c(-max(c(df1$PC1, df1$PC2)), max(c(df1$PC1, df1$PC2)))) +
#   # xlim(c(-max(c(df2$PC1, df2$PC2)), max(c(df2$PC1, df2$PC2)))) + ylim(c(-max(c(df2$PC1, df2$PC2)), max(c(df2$PC1, df2$PC2)))) +
#   theme_bw() +
#   # scale_alpha_manual(name = "Dataset", labels = factor(c("Training", "Validation"), levels = c("Training", "Validation")), guide = "legend", values = c(0.2, 0.8)) +
#   scale_color_manual(name = "Dataset", labels = factor(c("Training", "Validation"), levels = c("Training", "Validation")), guide = "legend", values = c("black", brewer.pal(n = 8, name = "Dark2")[6])) +
#   theme(panel.grid = element_blank()) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) 
# 
# 
# df1 <- merge(df1, combined.df.select, by = "Date.Time")
# summary(lm(df1$aou~df1$PC1))
# summary(lm(df1$temperature.sccoos~df1$PC1))
# summary(lm(df1$salinity~df1$PC1))
# summary(lm(df1$temperature.sccoos~df1$PC2))
# summary(lm(df1$salinity~df1$PC2))
# summary(lm(df1$DPD~df1$PC2))
# 
# a <- ggplot(df1, aes(x=PC1, y=PC2)) + 
#   geom_point(aes(color = aou)) +
#   geom_hline(yintercept=0, linetype="dotted") +
#   geom_vline(xintercept=0, linetype="dotted") +
#   coord_fixed() +
#   geom_segment(data=df2, aes(x=0, xend=PC1, y=0, yend=PC2),
#                color="black", arrow=arrow(length=unit(0.01,"npc"))) +
#   geom_text(data=df2,
#             aes(x=PC1,y=PC2,label=rownames(df2)),
#             # hjust= (1*sign(PC1)), vjust= (-0.5*sign(PC1))),
#             color="black", size=4) +
#   # scale_color_manual(values = c("grey69", "blue")) +
#   # scale_alpha_manual(values = c(0.1, 0.8)) +
#   scale_color_viridis_c() +
#   xlim(c(-max(c(df2$PC1, df2$PC2)), max(c(df2$PC1, df2$PC2)))) + ylim(c(-max(c(df2$PC1, df2$PC2)), max(c(df2$PC1, df2$PC2)))) +
#   theme_bw()
# 
# ggplotly(a)
# 

# ---- train partial model ----

set.seed(1234)

combined.gbm <- gbm(BOU ~ .,
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
  labs(x = "GBM Predictor", y = "Relative Influence") +
  scale_x_discrete(labels = rev(as.character(combined.gbm.summary$var))) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12))

sum(combined.gbm.summary$rel.inf[1:3])

## assess model fit ##
plot(combined.gbm$fit ~ combined.df.train$BOU,
     ylab = 'Predicted',
     xlab = 'Observed')
abline(0, 1, col = 'red')

summary(lm(combined.gbm$fit ~ combined.df.train$BOU))

plot((combined.gbm$fit - combined.df.train$BOU),
     type = 'h')



# ---- apply model to validation data ----

set.seed(1234)
BOU.est <- predict(combined.gbm, combined.df.test)

BOU.est.delta <- BOU.est - combined.df.test$BOU #difference between predicted and o2bio
AOU.BOU.delta <- combined.df.test$aou - combined.df.test$BOU #difference between uncorrected aou and o2bio

m1.RMSE.test <- sqrt(mean((BOU.est - combined.df.test$BOU)^2))
m1.perc.err.test <- median(abs(BOU.est - combined.df.test$BOU)/abs(combined.df.test$BOU))
m1.R2.test <- summary(lm(BOU.est ~ combined.df.test$BOU))$r.squared

BOU.est.full <- predict(combined.gbm, combined.df.select)
m1.RMSE.full <- sqrt(mean((BOU.est.full - combined.df.select$BOU)^2))
m1.perc.err.full <- mean(abs(BOU.est.full - combined.df.select$BOU)/abs(combined.df.select$BOU))
m1.R2.full <- summary(lm(BOU.est.full ~ combined.df.select$BOU))$r.squared

# # # ---- PCA of validation data ----
# # 
# # 
# # my.date.labels <- combined.df.test$Date.Time
# # env.pred.matrix <- as.matrix(combined.df.test[,predictors])
# # env.pred.matrix <- env.pred.matrix[,-9]
# # 
# # PCA <- rda(env.pred.matrix, scale = FALSE)
# # 
# # barplot(as.vector(PCA$CA$eig)/sum(PCA$CA$eig))
# # sum((as.vector(PCA$CA$eig)/sum(PCA$CA$eig))[1:2])
# # smry <- summary(PCA)
# # df1  <- data.frame(smry$sites[,1:2])# PC1 and PC2
# # df1$Date.Time <- my.date.labels
# # df2  <- data.frame(smry$species[,1:2])     # loadings for PC1 and PC2
# # 
# # # df2 <- df2[which((abs(df2$PC1) + abs(df2$PC2)) >= 0.3),]
# # 
# # df1 <- cbind(df1, combined.df.test)
# # df1$BOU.est <- BOU.est
# # df1$delta <- df1$BOU.est - df1$aou
# # df1$error <- df1$BOU.est - df1$BOU
# # 
# # ggplot(df1, aes(x=PC1, y=PC2)) + 
# #   geom_point(aes(color = error)) +
# #   geom_hline(yintercept=0, linetype="dotted") +
# #   geom_vline(xintercept=0, linetype="dotted") +
# #   coord_fixed() +
# #   # geom_segment(data=df2, aes(x=0, xend=PC1, y=0, yend=PC2),
# #   #              color="black", arrow=arrow(length=unit(0.01,"npc"))) +
# #   # geom_text(data=df2,
# #   #           aes(x=PC1,y=PC2,label=rownames(df2)),
# #   #               # hjust= (1*sign(PC1)), vjust= (-0.5*sign(PC1))),
# #   #           color="black", size=4) +
# #   xlim(c(-max(c(df1$PC1, df1$PC2)), max(c(df1$PC1, df1$PC2)))) + ylim(c(-max(c(df1$PC1, df1$PC2)), max(c(df1$PC1, df1$PC2)))) +
# #   # xlim(c(-max(c(df2$PC1, df2$PC2)), max(c(df2$PC1, df2$PC2)))) + ylim(c(-max(c(df2$PC1, df2$PC2)), max(c(df2$PC1, df2$PC2)))) +
# #   theme_bw() +
# #   scale_alpha_manual(name = "Dataset", labels = factor(c("Training", "Validation"), levels = c("Training", "Validation")), guide = "legend", values = c(0.2, 0.8)) +
# #   # scale_color_manual(name = "Dataset", labels = factor(c("Training", "Validation"), levels = c("Training", "Validation")), guide = "legend", values = c("grey69", "black")) +
# #   scale_color_viridis_c() +
# #   # theme(panel.grid = element_blank()) +
# #   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
# #   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) 
# # 
# # 
# # summary(lm(df1$delta~df1$PC1))
# # summary(lm(df1$error~df1$PC1))
# # summary(lm(df1$temperature.sccoos~df1$PC1))
# # summary(lm(df1$salinity~df1$PC1))
# # summary(lm(df1$temperature.sccoos~df1$PC2))
# # summary(lm(df1$salinity~df1$PC2))
# # 
# # 
# # a <- ggplot(df1, aes(x=PC1, y=PC2)) + 
# #   geom_point(aes(color = aou)) +
# #   geom_hline(yintercept=0, linetype="dotted") +
# #   geom_vline(xintercept=0, linetype="dotted") +
# #   coord_fixed() +
# #   geom_segment(data=df2, aes(x=0, xend=PC1, y=0, yend=PC2),
# #                color="black", arrow=arrow(length=unit(0.01,"npc"))) +
# #   geom_text(data=df2,
# #             aes(x=PC1,y=PC2,label=rownames(df2)),
# #             # hjust= (1*sign(PC1)), vjust= (-0.5*sign(PC1))),
# #             color="black", size=4) +
# #   # scale_color_manual(values = c("grey69", "blue")) +
# #   # scale_alpha_manual(values = c(0.1, 0.8)) +
# #   scale_color_viridis_c() +
# #   xlim(c(-max(c(df2$PC1, df2$PC2)), max(c(df2$PC1, df2$PC2)))) + ylim(c(-max(c(df2$PC1, df2$PC2)), max(c(df2$PC1, df2$PC2)))) +
# #   theme_bw()
# # 
# # ggplotly(a)
# # 
# # 
# # 
# 
# # ---- parameter hypertuning ----
# 
# hyper.grid <- expand.grid(
#   interaction.depth = 2:6,
#   n.trees = seq(1000, 4000, 1000),
#   bag.fraction = seq(0.3, 0.7, 0.1),
#   shrinkage = c(0.001, 0.01, 0.1),
#   RMSE = NA
# )
# 
# set.seed(1234)
# 
# for(i in 1:nrow(hyper.grid)){
#   
#   ## try clause necessary because some parameter combinations are
#   ## incompatible
#   
#   try({
#     
#     temp.model <- gbm(
#       formula = BOU ~ .,
#       data = combined.df.train,
#       interaction.depth = hyper.grid$interaction.depth[i],
#       distribution = 'gaussian',
#       n.trees       = hyper.grid$n.trees[i],
#       bag.fraction  = hyper.grid$bag.fraction[i],
#       shrinkage   = hyper.grid$shrinkage[i],
#       cv.folds = 2
#     )
#     
#     temp.predict <- predict(temp.model, combined.df.test)
#     temp.delta <- temp.predict - combined.df.test$BOU
#     hyper.grid$RMSE[i] <- sqrt(mean(temp.delta^2))
#   }, silent = F)
#   
#   print(paste(i, 'out of', nrow(hyper.grid), hyper.grid$RMSE[i]))
# }
# 
# hyper.grid <- na.omit(hyper.grid)
# 
# hist(hyper.grid$RMSE, breaks = 100)
# 
# final.gbm <- gbm(
#   formula = BOU ~ .,
#   data = combined.df.train,
#   interaction.depth = hyper.grid$interaction.depth[which.min(hyper.grid$RMSE)],
#   distribution = 'gaussian',
#   n.trees       = hyper.grid$n.trees[which.min(hyper.grid$RMSE)],
#   bag.fraction  = hyper.grid$bag.fraction[which.min(hyper.grid$RMSE)],
#   shrinkage   = hyper.grid$shrinkage[which.min(hyper.grid$RMSE)],
#   cv.folds = 2
# )
# 
# saveRDS(final.gbm, file = "2024-08-26_final_gbm.rds")
# 
# # ---- apply final model to test data ----
# 
# final.gbm <- readRDS("2024-08-26_final_gbm.rds")
# 
# BOU.est.final <- predict(final.gbm, combined.df.test)
# 
# plot.2b <- ggplot() +
#   geom_hline(aes(yintercept = 0), alpha = 0.2) +
#   geom_line(aes(x = combined.df.test.Date.Time, y = combined.df.test$BOU, color = "BOU"), lwd = 1, alpha = 0.7) +
#   geom_line(aes(x = combined.df.test.Date.Time, y = combined.df.test$aou, color = "AOU"), lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df.test.Date.Time, y = BOU.est.final, color = "Estimated BOU"), lwd = 1, alpha = 0.7) +
#   labs(x = "Date (2021)", y = "Oxygen Anomaly [μM]", color = "") +
#   ylim(c(-100, 100)) +
#   theme_bw() +
#   # scale_color_manual(name = "", labels = factor(c("AOU", "BOU", "Estimated BOU"), levels = c("BOU", "AOU", "Estimated BOU")), guide = "legend", values = c("BOU" = col2.o2bio, "AOU" = col1.aou, "Estimated BOU" = col5.other)) +
#   scale_color_brewer(palette = "Dark2") +
#   scale_x_datetime(date_breaks = "3 days", date_labels = "%b %d") +
#   theme(panel.grid = element_blank()) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   # theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   theme(legend.position = "none") 
# 
# plot.2c <- ggplot() +
#   geom_hline(aes(yintercept = 0), alpha = 0.2) +
#   geom_line(aes(x = combined.df.other.Date.Time, y = combined.df.other$BOU, color = "BOU"), lwd = 1, alpha = 0.7) +
#   geom_line(aes(x = combined.df.other.Date.Time, y = combined.df.other$aou, color = "AOU"), lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df.other.Date.Time, y = BOU.est.final, color = "Estimated BOU"), lwd = 1, alpha = 0.7) +
#   labs(x = "Date (2022)", y = "Oxygen Anomaly [μM]", color = "") +
#   # ylim(c(-100, 100)) +
#   theme_bw() +
#   # scale_color_manual(name = "", labels = factor(c("AOU", "BOU", "Estimated BOU"), levels = c("BOU", "AOU", "Estimated BOU")), guide = "legend", values = c("BOU" = col2.o2bio, "AOU" = col1.aou, "Estimated BOU" = col5.other)) +
#   scale_color_brewer(palette = "Dark2") +
#   scale_x_datetime(date_breaks = "3 days", date_labels = "%b %d") +
#   theme(panel.grid = element_blank()) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   # theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   theme(legend.position = "none")
# 
# 
# 
# plot.3b <- ggplot() +
#   geom_hline(aes(yintercept = 0), alpha = 0.2) +
#   geom_line(aes(x = combined.df.test.Date.Time, y = combined.df.test$BOU, color = "BOU"), lwd = 1, alpha = 0.7) +
#   geom_line(aes(x = combined.df.test.Date.Time, y = combined.df.test$aou, color = "AOU"), lwd = 1, alpha = 0.7) +
#   geom_line(aes(x = combined.df.test.Date.Time, y = BOU.est.final, color = "Estimated BOU"), lwd = 1, alpha = 0.7) +
#   labs(x = "Date", y = "Oxygen Anomaly [μM]", color = "") +
#   ylim(c(-100, 100)) +
#   theme_bw() +
#   # scale_color_manual(name = "", labels = factor(c("AOU", "BOU", "Estimated BOU"), levels = c("BOU", "AOU", "Estimated BOU")), guide = "legend", values = c("BOU" = col2.o2bio, "AOU" = col1.aou, "Estimated BOU" = col5.other)) +
#   scale_color_brewer(palette = "Dark2") +
#   scale_x_datetime(date_breaks = "3 days", date_labels = "%b %d") +
#   theme(panel.grid = element_blank()) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) 
# 
# # ---- calculate model performance on validation data ----
# 
# m2.RMSE.test <- sqrt(mean((BOU.est.final - combined.df.test$BOU)^2))
# m2.perc.err.test <- median(abs(BOU.est.final - combined.df.test$BOU)/abs(combined.df.test$BOU))
# m2.R2.test <- summary(lm(BOU.est.final ~ combined.df.test$BOU))$r.squared
# 
# BOU.est.final.full <- predict(final.gbm, combined.df.select)
# m2.RMSE.full <- sqrt(mean((BOU.est.final.full - combined.df.select$BOU)^2))
# m2.perc.err.full <- mean(abs(BOU.est.final.full - combined.df.select$BOU)/abs(combined.df.select$BOU))
# m2.R2.full <- summary(lm(BOU.est.final.full ~ combined.df.select$BOU))$r.squared
# 
# 
# ## various analyses of error reduction from the model
# # summary(model)
# # summary(model.old)
# 
# t.test(x = (BOU.est.final-combined.df.test$BOU), mu = 0)
# t.test(x = (combined.df.test$aou-combined.df.test$BOU), mu = 0)
# 
# hist(BOU.est.final-combined.df.test$BOU)
# hist(combined.df.test$aou-combined.df.test$BOU)
# 
# sum(abs(BOU.est.final-combined.df.test$BOU))
# sum(abs(combined.df.test$aou-combined.df.test$BOU))
# 
# plot.3c <- ggplot() +
#   geom_point(aes(x = combined.df.test$BOU, y = combined.df.test$aou, color = "AOU"), size = 2) +
#   geom_point(aes(combined.df.test$BOU, y = BOU.est.final, color = "Estimated BOU"), size = 2) +
#   geom_abline(slope = 1, intercept = 0, linewidth = 2, col= brewer.pal(n = 8, name = "Dark2")[2]) +
#   scale_y_continuous(breaks = c(-75, -50, -25, 0, 25)) +
#   labs(x = "Measured BOU [μM]", y = "Oxygen Anomaly [μM]", color = "") +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   # scale_color_brewer(palette = "Dark2") +
#   scale_color_manual(name = "", labels = factor(c("AOU", "Estimated BOU"), levels = c("AOU", "Estimated BOU")), guide = "legend", values = c("AOU" = brewer.pal(n = 8, name = "Dark2")[1], "Estimated BOU" = brewer.pal(n = 8, name = "Dark2")[3])) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) 
# 
# png("../Figures/manuscript_figures/2024-08-26_fig2.png", width = 700, height = 500)
# 
# (plot.2a) / (plot.2b + plot.2c)
# 
# dev.off()
# 
# png("../Figures/manuscript_figures/2024-08-26_fig3.png", width = 1000, height = 250)
# 
# plot.3b + plot.3c + plot_layout(widths = c(2,1))
# 
# dev.off()
# 
# # ---- full model ----
# 
# set.seed(1234)
# 
# full.gbm <- gbm(BOU ~ .,
#                 data = na.omit(combined.df.select[,predictors]),
#                 interaction.depth = hyper.grid$interaction.depth[which.min(hyper.grid$RMSE)],
#                 distribution = 'gaussian',
#                 n.trees       = hyper.grid$n.trees[which.min(hyper.grid$RMSE)],
#                 bag.fraction  = hyper.grid$bag.fraction[which.min(hyper.grid$RMSE)],
#                 shrinkage   = hyper.grid$shrinkage[which.min(hyper.grid$RMSE)],
#                 cv.folds = 2)
# 
# 
# set.seed(1234)
# BOU.est.final <- predict(full.gbm, combined.df.select)
# 
# m3.RMSE.test <- sqrt(mean((BOU.est.final - combined.df.select$BOU)^2))
# m3.perc.err.test <- median(abs(BOU.est.final - combined.df.select$BOU)/abs(combined.df.select$BOU))
# m3.R3.test <- summary(lm(BOU.est.final ~ combined.df.select$BOU))$r.squared
# 
# 
# save(list = c('combined.gbm', 'full.gbm', 'hyper.grid'), file = '20240826_model.Rdata')
# 
# 
# 
# # ---- apply model to full miniDot timeseries ----
# 
# combined.df <- readRDS("2024-08-26_combined_env_data_hourly.rds")
# load(file = "20240826_model.Rdata")
# 
# combined.df$BOU.estimated <- predict(full.gbm, combined.df[,c(full.gbm$var.names, 'Date.Time')])
# 
# combined.df.select <- na.omit(combined.df)
# mean(abs(combined.df.select$BOU - combined.df.select$BOU.estimated)/abs(combined.df.select$BOU)) # perc error, mean
# median(abs(combined.df.select$BOU - combined.df.select$BOU.estimated)/abs(combined.df.select$BOU)) # perc error, median
# 
# sqrt(mean((combined.df.select$BOU - combined.df.select$BOU.estimated)^2)) #RMSE
# # mean(abs(combined.df.select$BOU - combined.df.select$BOU.estimated)/abs(combined.df.select$BOU))
# 
# mean(combined.df.select$BOU.estimated)
# sd(combined.df.select$BOU.estimated)
# 
# 
# 
# 
# 
# ## assess variable significance ##
# full.gbm.summary <- summary(full.gbm) 
# 
# full.gbm.summary$var <- factor(full.gbm.summary$var, levels = rev(rownames(full.gbm.summary)))
# ggplot(data = full.gbm.summary) +
#   geom_col(aes(x = var, y = rel.inf), fill = c(col5.other, col5.other, col5.other, col1.aou, col1.aou, col1.aou, col1.aou, col1.aou)) +
#   labs(x = "GBM Predictor", y = "Relative Influence") +
#   scale_x_discrete(labels = rev(as.character(full.gbm.summary$var))) +
#   coord_flip() +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12))
# 
# sum(full.gbm.summary$rel.inf[1:3])
# 
# saveRDS(combined.df, "2024-08-09_bou_est_df.rds")
# 
# # ---- aggregate by day if desired ----
# 
# combined.df$Date <- parse_date_time(paste(day(combined.df$Date.Time), month(combined.df$Date.Time), year(combined.df$Date.Time), sep = "-"), orders = "dmY")
# my.colnames <- colnames(combined.df)[-which(colnames(combined.df) == "Date.Time")]
# combined.df.daily <- combined.df %>% group_by(Date) %>% summarize_all(mean)
# saveRDS(combined.df.daily, file = "2024-08-09_combined_BOUest_env_daily.rds")
# 
# # 
# # ggplot(data = combined.df.daily) +
# #   geom_hline(aes(yintercept = 0), alpha = 0.5) +
# #   geom_line(aes(x = Date.Time, y = aou, color = "AOU", group = plot.group), lwd = 1, alpha = 0.7) +
# #   geom_line(aes(x = Date.Time, y = BOU, color = "BOU", group = plot.group), lwd = 1, alpha = 0.7) +
# #   geom_line(aes(x = Date.Time, y = BOU.estimated, color = "Estimated BOU", group = plot.group), lwd = 1, alpha = 0.7) +
# #   geom_hline(aes(yintercept = mean(aou), color = "AOU")) +
# #   geom_hline(aes(yintercept = mean(na.omit(BOU)), color = "BOU")) +
# #   geom_hline(aes(yintercept = mean(BOU.estimated), color = "Estimated BOU")) +
# #   labs(x = "Date", y = "Oxygen Anomaly [μM]", color = "") +
# #   theme_bw() +
# #   ylim(c(-250, 250)) +
# #   scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
# #   # theme(panel.grid = element_blank()) +
# #   scale_color_brewer(palette = "Dark2") +
# #   # scale_color_manual(name = "", labels = factor(c("AOU", "BOU", "Estimated BOU"), levels = c("AOU", "BOU", "Estimated BOU")), guide = "legend", values = c("BOU" = col2.o2bio, "AOU" = col1.aou, "Estimated BOU" = col5.other)) +
# #   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
# #   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
# #   theme(panel.grid = element_blank())
# 
# # ---- analysis of temporal trends (and visualization) ----
# 
# combined.df.select
# 
# ## summary of uncorrected AOU
# summary(combined.df$aou) # full time series, not just what is used in the model
# sd(combined.df$aou)
# t.test(combined.df$aou, mu = 0)
# 
# ## summary of estimated BOU
# summary(combined.df$BOU.estimated) # full time series, not just what is used in the model
# sd(combined.df$BOU.estimated)
# t.test(combined.df$BOU.estimated, mu = 0)
# 
# ## summary of measured BOU
# summary(combined.df$BOU) # full time series, not just what is used in the model
# sd(na.omit(combined.df$BOU))
# t.test(combined.df$BOU, mu = 0)
# 
# 
# combined.df <- set.plot.groups(combined.df)
# 
# ## full corrected aou time series (with uncorrected aou and [o2]bio)
# ggplot(data = combined.df) +
#   geom_hline(aes(yintercept = 0), alpha = 0.5) +
#   geom_line(aes(x = Date.Time, y = aou, color = "AOU", group = plot.group), lwd = 1, alpha = 0.7) +
#   geom_line(aes(x = Date.Time, y = BOU, color = "BOU", group = plot.group), lwd = 1, alpha = 0.7) +
#   geom_line(aes(x = Date.Time, y = BOU.estimated, color = "Estimated BOU", group = plot.group), lwd = 1, alpha = 0.7) +
#   labs(x = "Date", y = "Oxygen Anomaly [μM]", color = "") +
#   theme_bw() +
#   ylim(c(-250, 250)) +
#   scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
#   # theme(panel.grid = element_blank()) +
#   scale_color_brewer(palette = "Dark2") +
#   # scale_color_manual(name = "", labels = factor(c("AOU", "BOU", "Estimated BOU"), levels = c("AOU", "BOU", "Estimated BOU")), guide = "legend", values = c("BOU" = col2.o2bio, "AOU" = col1.aou, "Estimated BOU" = col5.other)) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   theme(panel.grid = element_blank())
# 
# 
# mean(combined.df$aou)
# sd(combined.df$aou)
# t.test(combined.df$aou, mu = 0)
# 
# mean(na.omit(combined.df$BOU))
# sd(na.omit(combined.df$BOU))
# t.test(na.omit(combined.df$BOU), mu = 0)
# 
# mean(na.omit(combined.df$BOU.estimated))
# sd(na.omit(combined.df$BOU.estimated))
# t.test(na.omit(combined.df$BOU.estimated), mu = 0)
# 
# plot.5a <- ggplot(data = combined.df) +
#   geom_hline(aes(yintercept = 0), alpha = 0.5) +
#   geom_line(aes(x = Date.Time, y = aou, color = "AOU", group = plot.group), lwd = 1, alpha = 0.7) +
#   geom_line(aes(x = Date.Time, y = BOU, color = "BOU", group = plot.group), lwd = 1, alpha = 0.7) +
#   geom_line(aes(x = Date.Time, y = BOU.estimated, color = "Estimated BOU", group = plot.group), lwd = 1, alpha = 0.7) +
#   geom_hline(aes(yintercept = mean(aou), color = "AOU"), linewidth = 1) +
#   geom_hline(aes(yintercept = mean(na.omit(BOU)), color = "BOU"), linewidth = 1) +
#   geom_hline(aes(yintercept = mean(BOU.estimated), color = "Estimated BOU"), linewidth = 1) +
#   labs(x = "Date", y = "Oxygen Anomaly [μM]", color = "") +
#   theme_bw() +
#   ylim(c(-250, 250)) +
#   scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
#   # theme(panel.grid = element_blank()) +
#   scale_color_brewer(palette = "Dark2") +
#   # scale_color_manual(name = "", labels = factor(c("AOU", "BOU", "Estimated BOU"), levels = c("AOU", "BOU", "Estimated BOU")), guide = "legend", values = c("BOU" = col2.o2bio, "AOU" = col1.aou, "Estimated BOU" = col5.other)) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   theme(panel.grid = element_blank())
# 
# test.data <- combined.df[which(combined.df$Date.Time %in% combined.df.test.Date.Time),]
# plot.5b <- ggplot(data = test.data) +
#   geom_hline(aes(yintercept = 0), alpha = 0.2) +
#   geom_line(aes(x = Date.Time, y = BOU, color = "BOU"), lwd = 1, alpha = 0.7) +
#   geom_line(aes(x = Date.Time, y = aou, color = "AOU"), lwd = 1, alpha = 0.7) +
#   geom_line(aes(x = Date.Time, y = BOU.estimated, color = "Estimated BOU"), lwd = 1, alpha = 0.7) +
#   labs(x = "Date (2021)", y = "Oxygen Anomaly [μM]", color = "") +
#   ylim(c(-100, 100)) +
#   theme_bw() +
#   # scale_color_manual(name = "", labels = factor(c("AOU", "BOU", "Estimated BOU"), levels = c("BOU", "AOU", "Estimated BOU")), guide = "legend", values = c("BOU" = col2.o2bio, "AOU" = col1.aou, "Estimated BOU" = col5.other)) +
#   scale_color_brewer(palette = "Dark2") +
#   scale_x_datetime(date_breaks = "3 days", date_labels = "%b %d") +
#   theme(panel.grid = element_blank()) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.position = "none") +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) 
# 
# other.data <- combined.df[which(combined.df$Date.Time %in% combined.df.other.Date.Time),]
# plot.5c <- ggplot(other.data) +
#   geom_hline(aes(yintercept = 0), alpha = 0.2) +
#   geom_line(aes(x = Date.Time, y = BOU, color = "BOU"), lwd = 1, alpha = 0.7) +
#   geom_line(aes(x = Date.Time, y = aou, color = "AOU"), lwd = 1, alpha = 0.7) +
#   geom_line(aes(x = Date.Time, y = BOU.estimated, color = "Estimated BOU"), lwd = 1, alpha = 0.7) +
#   labs(x = "Date", y = "Oxygen Anomaly [μM]", color = "") +
#   # ylim(c(-100, 100)) +
#   theme_bw() +
#   # scale_color_manual(name = "", labels = factor(c("AOU", "BOU", "Estimated BOU"), levels = c("BOU", "AOU", "Estimated BOU")), guide = "legend", values = c("BOU" = col2.o2bio, "AOU" = col1.aou, "Estimated BOU" = col5.other)) +
#   scale_color_brewer(palette = "Dark2") +
#   scale_x_datetime(date_breaks = "3 days", date_labels = "%b %d") +
#   theme(panel.grid = element_blank()) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.position = "none") +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) 
# 
# png("../Figures/manuscript_figures/2024-08-09_fig5.png", width = 1000, height = 500)
# 
# plot.5a / (plot.5b + plot.5c)
# 
# dev.off()
# 
# 
# ## ---- calculate differences between uaou, caou, and [o2]bio ----
# combined.df$delta <- (combined.df$aou - combined.df$BOU.estimated)/combined.df$BOU
# combined.df$error.c <- (combined.df$BOU.estimated - combined.df$BOU)/combined.df$BOU
# combined.df$error.c.ab <- combined.df$BOU.estimated - combined.df$BOU
# combined.df$error.u <- (combined.df$aou - combined.df$BOU)/combined.df$BOU
# 
# combined.df$Year <- year(combined.df$Date.Time)
# combined.df$date.mm.dd <- parse_date_time(paste(month(combined.df$Date.Time), day(combined.df$Date.Time), hour(combined.df$Date.Time), minute(combined.df$Date.Time), second(combined.df$Date.Time),sep = "-"), orders = "mdHMS")
# 
# summary(abs(combined.df$error.u))
# summary(abs(combined.df$error.c))
# summary(lm(combined.df$delta~combined.df$error.u))
# 
# try.it <- combined.df[which(combined.df$Year >= 2021),]
# ggplot(data = try.it) +
#   geom_hline(aes(yintercept = 0), alpha = 0.5) +
#   # geom_line(aes(x = date.mm.dd, y = delta), color = col5.other, lwd = 1, alpha = 0.7) +
#   geom_line(aes(x = date.mm.dd, y = (error.c)), color = col3.BOU_est, lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df$date.mm.dd, y = abs(combined.df$error.c)), color = "blue", lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = date.mm.dd, y = (error.u)), color = col1.aou, lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df$date.mm.dd, y = abs(combined.df$error.c)), color = "blue", lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df$date.mm.dd, y = (combined.df$aou)), color = "black", lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df$date.mm.dd, y = (combined.df$delta), color = "Delta"), lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df$date.mm.dd, y = (combined.df$delta2)), color = "red", lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df$date.mm.dd, y = (combined.df$delta3)), color = "blue", lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df$date.mm.dd, y = (combined.df$BOU)), color = "blue", lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df$date.mm.dd, y = (combined.df$temperature.sccoos*3)), color = "black", lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df.select$Date.Time, y = combined.df.select$BOU, color = "[O2]bio"), lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df$Date.Time, y = combined.df$BOU.estrected, color = "Corrected aou"), lwd = 1, alpha = 0.7) +
#   labs(x = "Date", y = "Percent Change between BOU and Estimated BOU [%]") +
#   theme_bw() +
#   # ylim(c(-250, 250)) +
#   # scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
#   # theme(panel.grid = element_blank()) +
#   scale_color_manual(name = "", labels = factor(c("Delta"), levels = c("Delta")), guide = "legend", values = c("Delta" = col5.other)) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   facet_wrap(.~Year, ncol = 1)
# 
# # plot(x = try.it$error.c, y = try.it$error.u)
# 
# sccoos.hourly$Date.Time <- parse_date_time(sccoos.hourly$Date.Time, orders = "Ymd HMS")
# try.it <- merge(x = combined.df, y = sccoos.hourly[,c(1,5)], by = "Date.Time", all.x = T, all.y = F)
# 
# ggplot(data = try.it) +
#   # geom_line(aes(x = date.mm.dd, y = chlorophyll), color = "blue") +
#   # geom_line(aes(x = date.mm.dd, y = log(abs(error.c))), color = "red") +
#   theme_bw() +
#   # facet_wrap(.~Year, ncol = 1)
#   geom_point(aes(x = abs(error.c), y = chlorophyll), color = "darkgreen")
# 
# summary(lm(formula = try.it$error.c ~ try.it$chlorophyll))
# 
# 
# summary(lm(formula = combined.df$delta ~ combined.df$aou))
# plot(combined.df$delta ~ combined.df$aou)
# summary(lm(formula = combined.df$error.u ~ combined.df$aou))
# summary(lm(formula = combined.df$error.c ~ combined.df$aou))
# 
# summary(lm(formula = combined.df$delta ~ combined.df$BOU.estrected))
# summary(lm(formula = combined.df$delta ~ combined.df$BOU))
# 
# ggplot(data = combined.df) +
#   geom_hline(aes(yintercept = 0), alpha = 0.5) +
#   geom_line(aes(x = combined.df$date.mm.dd, y = combined.df$temperature.sccoos, color = "Delta"), lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df.select$Date.Time, y = combined.df.select$BOU, color = "[O2]bio"), lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df$Date.Time, y = combined.df$BOU.estrected, color = "Corrected aou"), lwd = 1, alpha = 0.7) +
#   labs(x = "Date", y = "Delta [μM]") +
#   theme_bw() +
#   # ylim(c(-250, 250)) +
#   # scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
#   # theme(panel.grid = element_blank()) +
#   scale_color_manual(name = "", labels = factor(c("Delta"), levels = c("Delta")), guide = "legend", values = c("Delta" = col5.other)) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   facet_wrap(.~Year, ncol = 1)
# range(combined.df$temperature.sccoos)
# O2sat(33.5, min(combined.df$temperature.sccoos))
# O2sat(33.5, max(combined.df$temperature.sccoos))
# 
# 
# combined.df$season <- "winter"
# combined.df$season[which(combined.df$date.mm.dd > parse_date_time("04-01", orders = "md") & combined.df$date.mm.dd < parse_date_time("10-01", orders = "md"))] <- "summer"
# 
# ggplot(data = combined.df) +
#   geom_hline(aes(yintercept = 0), alpha = 0.5) +
#   geom_boxplot(aes(x = combined.df$season, y = abs(combined.df$error.c)), lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df.select$Date.Time, y = combined.df.select$BOU, color = "[O2]bio"), lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df$Date.Time, y = combined.df$BOU.estrected, color = "Corrected aou"), lwd = 1, alpha = 0.7) +
#   labs(x = "Date", y = "Perc. Error [%]") +
#   theme_bw() +
#   # ylim(c(-250, 250)) +
#   # scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
#   # theme(panel.grid = element_blank()) +
#   # scale_color_manual(name = "", labels = factor(c("Delta"), levels = c("Delta")), guide = "legend", values = c("Delta" = col5.other)) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) 
# 
# t.test(x = abs(combined.df$error.c[which(combined.df$season == "summer")]), y = abs(combined.df$error.c[which(combined.df$season == "winter")]))
# 
# ggplot(data = combined.df) +
#   geom_hline(aes(yintercept = 0), alpha = 0.5) +
#   geom_boxplot(aes(x = combined.df$season, y = abs(combined.df$error.c.ab)), lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df.select$Date.Time, y = combined.df.select$BOU, color = "[O2]bio"), lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df$Date.Time, y = combined.df$BOU.estrected, color = "Corrected aou"), lwd = 1, alpha = 0.7) +
#   labs(x = "Date", y = "Perc. Error [%]") +
#   theme_bw() +
#   # ylim(c(-250, 250)) +
#   # scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
#   # theme(panel.grid = element_blank()) +
#   # scale_color_manual(name = "", labels = factor(c("Delta"), levels = c("Delta")), guide = "legend", values = c("Delta" = col5.other)) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) 
# 
# t.test(x = abs(combined.df$error.c.ab[which(combined.df$season == "summer")]), y = abs(combined.df$error.c.ab[which(combined.df$season == "winter")]))
# 
# 
# ggplot(data = combined.df) +
#   geom_hline(aes(yintercept = 0), alpha = 0.5) +
#   geom_boxplot(aes(x = combined.df$season, y = abs(combined.df$BOU)), lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df.select$Date.Time, y = combined.df.select$BOU, color = "[O2]bio"), lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df$Date.Time, y = combined.df$BOU.estrected, color = "Corrected aou"), lwd = 1, alpha = 0.7) +
#   labs(x = "Date", y = "BOU [uM]") +
#   theme_bw() +
#   # ylim(c(-250, 250)) +
#   # scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
#   # theme(panel.grid = element_blank()) +
#   # scale_color_manual(name = "", labels = factor(c("Delta"), levels = c("Delta")), guide = "legend", values = c("Delta" = col5.other)) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) 
# 
# t.test(x = abs(combined.df$BOU[which(combined.df$season == "summer")]), y = abs(combined.df$BOU[which(combined.df$season == "winter")]))
# 
# ggplot(data = combined.df) +
#   geom_hline(aes(yintercept = 0), alpha = 0.5) +
#   geom_boxplot(aes(x = combined.df$season, y = abs(combined.df$BOU.estimated)), lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df.select$Date.Time, y = combined.df.select$BOU, color = "[O2]bio"), lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df$Date.Time, y = combined.df$BOU.estrected, color = "Corrected aou"), lwd = 1, alpha = 0.7) +
#   labs(x = "Date", y = "Estimated BOU [uM]") +
#   theme_bw() +
#   # ylim(c(-250, 250)) +
#   # scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
#   # theme(panel.grid = element_blank()) +
#   # scale_color_manual(name = "", labels = factor(c("Delta"), levels = c("Delta")), guide = "legend", values = c("Delta" = col5.other)) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) 
# 
# t.test(x = abs(combined.df$BOU.estimated[which(combined.df$season == "summer")]), y = abs(combined.df$BOU.estimated[which(combined.df$season == "winter")]))
# 
# ggplot(data = combined.df) +
#   geom_hline(aes(yintercept = 0), alpha = 0.5) +
#   geom_boxplot(aes(x = combined.df$season, y = abs(combined.df$error.u)), lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df.select$Date.Time, y = combined.df.select$BOU, color = "[O2]bio"), lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df$Date.Time, y = combined.df$BOU.estrected, color = "Corrected aou"), lwd = 1, alpha = 0.7) +
#   labs(x = "Date", y = "Perc. Error [%]") +
#   theme_bw() +
#   # ylim(c(-250, 250)) +
#   # scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
#   # theme(panel.grid = element_blank()) +
#   # scale_color_manual(name = "", labels = factor(c("Delta"), levels = c("Delta")), guide = "legend", values = c("Delta" = col5.other)) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) 
# 
# t.test(x = abs(combined.df$error.u[which(combined.df$season == "summer")]), y = abs(combined.df$error.u[which(combined.df$season == "winter")]))
# 
# 
# ggplot(data = combined.df) +
#   geom_hline(aes(yintercept = 0), alpha = 0.5) +
#   geom_boxplot(aes(x = combined.df$season, y = abs(combined.df$error.c)), lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df.select$Date.Time, y = combined.df.select$BOU, color = "[O2]bio"), lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df$Date.Time, y = combined.df$BOU.estrected, color = "Corrected aou"), lwd = 1, alpha = 0.7) +
#   labs(x = "Date", y = "Delta [μM]") +
#   theme_bw() +
#   # ylim(c(-250, 250)) +
#   # scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
#   # theme(panel.grid = element_blank()) +
#   # scale_color_manual(name = "", labels = factor(c("Delta"), levels = c("Delta")), guide = "legend", values = c("Delta" = col5.other)) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) 
# 
# t.test(x = combined.df$error.c[which(combined.df$season == "summer")], y = combined.df$error.c[which(combined.df$season == "winter")])
# 
# ggplot(data = combined.df) +
#   geom_hline(aes(yintercept = 0), alpha = 0.5) +
#   geom_boxplot(aes(x = combined.df$season, y = abs(combined.df$error.u)), lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df.select$Date.Time, y = combined.df.select$BOU, color = "[O2]bio"), lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df$Date.Time, y = combined.df$BOU.estrected, color = "Corrected aou"), lwd = 1, alpha = 0.7) +
#   labs(x = "Date", y = "Delta [μM]") +
#   theme_bw() +
#   # ylim(c(-250, 250)) +
#   # scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
#   # theme(panel.grid = element_blank()) +
#   # scale_color_manual(name = "", labels = factor(c("Delta"), levels = c("Delta")), guide = "legend", values = c("Delta" = col5.other)) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) 
# 
# t.test(x = combined.df$error.u[which(combined.df$season == "summer")], y = combined.df$error.u[which(combined.df$season == "winter")])
# 
# 
# combined.df$delta <- abs(combined.df$BOU.estrected - combined.df$BOU)
# 
# ggplot(data = combined.df) +
#   geom_hline(aes(yintercept = 0), alpha = 0.5) +
#   geom_boxplot(aes(x = combined.df$season, y = abs(combined.df$delta)), lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df.select$Date.Time, y = combined.df.select$BOU, color = "[O2]bio"), lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df$Date.Time, y = combined.df$BOU.estrected, color = "Corrected aou"), lwd = 1, alpha = 0.7) +
#   labs(x = "Date", y = "Delta [μM]") +
#   theme_bw() +
#   # ylim(c(-250, 250)) +
#   # scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
#   # theme(panel.grid = element_blank()) +
#   # scale_color_manual(name = "", labels = factor(c("Delta"), levels = c("Delta")), guide = "legend", values = c("Delta" = col5.other)) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) 
# 
# t.test(x = combined.df$delta[which(combined.df$season == "summer")], y = combined.df$delta[which(combined.df$season == "winter")])
# 
# 
# 
# combined.df$delta <- abs(combined.df$aou - combined.df$BOU)
# 
# ggplot(data = combined.df) +
#   geom_hline(aes(yintercept = 0), alpha = 0.5) +
#   geom_boxplot(aes(x = combined.df$season, y = abs(combined.df$delta)), lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df.select$Date.Time, y = combined.df.select$BOU, color = "[O2]bio"), lwd = 1, alpha = 0.7) +
#   # geom_line(aes(x = combined.df$Date.Time, y = combined.df$BOU.estrected, color = "Corrected aou"), lwd = 1, alpha = 0.7) +
#   labs(x = "Date", y = "Delta [μM]") +
#   theme_bw() +
#   # ylim(c(-250, 250)) +
#   # scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
#   # theme(panel.grid = element_blank()) +
#   # scale_color_manual(name = "", labels = factor(c("Delta"), levels = c("Delta")), guide = "legend", values = c("Delta" = col5.other)) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) 
# 
# t.test(x = combined.df$delta[which(combined.df$season == "summer")], y = combined.df$delta[which(combined.df$season == "winter")])
# 
# 
# 
# 
# # temp <- merge(combined.df, combined.df.select, by = "Date.Time")
# # temp$delta <- temp$BOU - temp$BOU.estrected
# temp <- combined.df
# temp$delta <- temp$BOU.estrected - temp$aou
# temp$Year <- year(temp$Date.Time)
# temp$date.mm.dd <- parse_date_time(paste(month(temp$Date.Time), day(temp$Date.Time), sep = "-"), orders = "md")
# 
# ggplot() +
#   geom_line(aes(x = temp$Date.Time, y = temp$delta), lwd = 1, color = "blue") +
#   labs(x = "Date", y = "Delta [μM]") +
#   theme_bw() +
#   # ylim(c(-250, 250)) +
#   scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
#   # theme(panel.grid = element_blank()) +
#   geom_hline(yintercept = 0) +
#   # scale_color_manual(name = "", labels = factor(c("[O2]bio", "aou", "Corrected aou"), levels = c("[O2]bio", "aou", "Corrected aou")), guide = "legend", values = c("[O2]bio" = col2.o2bio, "aou" = col1.aou, "Corrected aou" = col3.BOU_est)) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   theme(panel.grid.major.x = element_line(color = "black", linewidth = 1))
# 
# ggplot(data = temp) +
#   geom_line(aes(x = date.mm.dd, y = delta), lwd = 1, color = "blue") +
#   labs(x = "Date", y = "Delta [μM]") +
#   theme_bw() +
#   # ylim(c(-250, 250)) +
#   # scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
#   # theme(panel.grid = element_blank()) +
#   geom_hline(yintercept = 0) +
#   # scale_color_manual(name = "", labels = factor(c("[O2]bio", "aou", "Corrected aou"), levels = c("[O2]bio", "aou", "Corrected aou")), guide = "legend", values = c("[O2]bio" = col2.o2bio, "aou" = col1.aou, "Corrected aou" = col3.BOU_est)) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   theme(panel.grid.major.x = element_line(color = "black", linewidth = 1)) +
#   facet_wrap(Year~., ncol = 1)
# 
# 
# # ---- GBM chunk cross validation ----
# 
# combined.df.select <- readRDS("2024-08-09_combined.df.select.rds")
# cross.val.df <- na.omit(combined.df.select[,c("Date.Time",predictors)])
# 
# # random.dates <- sample(cross.val.df$Date.Time[1:round(nrow(cross.val.df)*0.9)], size = 100, replace = F)
# 
# all.valid.dates <- c(cross.val.df$Date.Time[1])
# for(d in 1:nrow(cross.val.df)){
#   
#   date1 <- cross.val.df$Date.Time[d]
#   date2 <- date1 + 2*604800
#   
#   my.df <- cross.val.df$Date.Time[which(cross.val.df$Date.Time >= date1 & cross.val.df$Date.Time <= date2)]
#   
#   if(length(my.df) >= 0.8*336){ # number of hours in two weeks *0.9 in case a little missing data
#     
#     all.valid.dates <- c(all.valid.dates, date1)
#     
#   }
#   
# }
# 
# all.valid.dates <- all.valid.dates[-1]
# all.valid.days <- all.valid.dates[which(hour(all.valid.dates) == 0)]
# 
# my.cross.val.values <- as.data.frame(matrix(data = NA, nrow = length(all.valid.days), ncol = 15))
# colnames(my.cross.val.values) <- c("Date.Time", "n.days", "perc.data", "RMSE_no_model", "R2_no_model", "p_no_model", "avg_perc_error_no_model", "RMSE_model_hyper", "R2_model_hyper", "p_model_hyper", "avg_perc_error_hyper", "RMSE_model_hyper_full", "R2_model_hyper_full", "p_model_hyper_full", "avg_perc_error_hyper_full")
# 
# my.cross.val.values$Date.Time <- all.valid.days
# 
# for(d in 1:length(all.valid.days)){
#   
#   date1 <- all.valid.days[d]
#   date2 <- date1 + 2*604800
#   
#   index <- which(cross.val.df$Date.Time >= date1 & cross.val.df$Date.Time <= date2)
#   
#   cross.val.df.nd <- cross.val.df[,which(colnames(cross.val.df) != "Date.Time")]
#   
#   cross.val.test <- cross.val.df.nd[index,]
#   cross.val.train <- cross.val.df.nd[-index,] 
#   
#   my.cross.val.values$perc.data[d] <- nrow(cross.val.test)/nrow(cross.val.train)
#   
#   my.cross.val.values$RMSE_no_model[d] <- sqrt(mean((cross.val.test$BOU - cross.val.test$aou)^2))
#   my.cross.val.values$R2_no_model[d] <- summary(lm(formula = cross.val.test$BOU~cross.val.test$aou))$r.squared
#   my.cross.val.values$p_no_model[d] <- summary(lm(formula = cross.val.test$BOU~cross.val.test$aou))$coefficients[2,4]
#   my.cross.val.values$avg_perc_error_no_model[d] <- mean((cross.val.test$BOU - cross.val.test$aou)/cross.val.test$BOU)
#   
#   temp.gbm <- gbm(BOU ~ .,
#                   data = cross.val.train,
#                   interaction.depth = hyper.grid$interaction.depth[which.min(hyper.grid$RMSE)],
#                   distribution = 'gaussian',
#                   n.trees       = hyper.grid$n.trees[which.min(hyper.grid$RMSE)],
#                   bag.fraction  = hyper.grid$bag.fraction[which.min(hyper.grid$RMSE)],
#                   shrinkage   = hyper.grid$shrinkage[which.min(hyper.grid$RMSE)],
#                   cv.folds = 2)
#   
#   
#   temp.prediction <- predict(temp.gbm, cross.val.test)
#   
#   my.cross.val.values$RMSE_model_hyper[d] <- sqrt(mean((temp.prediction-cross.val.test$BOU)^2))
#   my.cross.val.values$R2_model_hyper[d] <- summary(lm(formula = temp.prediction~cross.val.test$BOU))$r.squared
#   my.cross.val.values$p_model_hyper[d] <- summary(lm(formula = temp.prediction~cross.val.test$BOU))$coefficients[2,4]
#   my.cross.val.values$avg_perc_error_hyper[d] <- mean((cross.val.test$BOU - temp.prediction)/cross.val.test$BOU)
#   
#   
#   full.prediction <- predict(temp.gbm, cross.val.df)
#   
#   my.cross.val.values$RMSE_model_hyper_full[d] <- sqrt(mean((full.prediction-cross.val.df$BOU)^2))
#   my.cross.val.values$R2_model_hyper_full[d] <- summary(lm(formula = full.prediction~cross.val.df$BOU))$r.squared
#   my.cross.val.values$p_model_hyper_full[d] <- summary(lm(formula = full.prediction~cross.val.df$BOU))$coefficients[2,4]
#   my.cross.val.values$avg_perc_error_hyper_full[d] <- mean((cross.val.df$BOU - full.prediction)/cross.val.df$BOU)
#   
#   
#   print(paste("Cross-Validation", d, "of", length(all.valid.days), "complete"))
#   
# }
# 
# saveRDS(my.cross.val.values, file = "2024-08-09_cross_validation_results.rds")
# my.cross.val.values <- readRDS("2024-08-09_cross_validation_results.rds")
# 
# hist(my.cross.val.values$RMSE_no_model)
# hist(my.cross.val.values$RMSE_model_hyper)
# hist(my.cross.val.values$RMSE_model_hyper_full)
# 
# hist(my.cross.val.values$avg_perc_error_no_model)
# hist(my.cross.val.values$avg_perc_error_hyper)
# hist(my.cross.val.values$avg_perc_error_hyper_full)
# 
# hist(my.cross.val.values$R2_no_model)
# hist(my.cross.val.values$R2_model_hyper)
# hist(my.cross.val.values$R2_model_hyper_full)
# 
# colnames(my.cross.val.values)[c(7,11,15)] <- c("No Model", "2-week Validation", "Full Time Series")
# # my.cross.val.values$`Full Time Series` <- 10*my.cross.val.values$`Full Time Series`
# # cross.val.long <- my.cross.val.values[,c(1,7,11,15)] %>% pivot_longer(cols = colnames(my.cross.val.values)[c(7,11,15)], names_to = "model", values_to = "my.perc.error")
# cross.val.long <- my.cross.val.values[,c(1,11)] %>% pivot_longer(cols = colnames(my.cross.val.values)[c(11)], names_to = "model", values_to = "my.perc.error")
# cross.val.long$Date <- parse_date_time(paste(month(cross.val.long$Date.Time), day(cross.val.long$Date.Time), sep = "-"), orders = "md")
# cross.val.long$Year <- year(cross.val.long$Date.Time)
# 
# ggplot(data = cross.val.long) +
#   # geom_line(aes(x = Date.Time, y = my.perc.error, color = model), size = 0.5) +
#   # geom_point(aes(x = Date.Time, y = my.perc.error, color = model), size = 2) +
#   geom_line(aes(x = Date.Time, y = abs(my.perc.error), color = model), size = 0.5) +
#   geom_point(aes(x = Date.Time, y = abs(my.perc.error), color = model), size = 2) +
#   scale_x_datetime(date_breaks = "1 month") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90))
# 
# 
# ggplot(data = cross.val.long) +
#   # geom_line(aes(x = Date, y = my.RMSE, color = model), size = 1) +
#   geom_point(aes(x = Date, y = abs(my.perc.error), color = factor(model, levels = c("No Model", "2-week Validation", "Full Time Series"))), size = 2) +
#   # geom_point(aes(x = Date, y = my.perc.error, color = model), size = 2) +
#   facet_wrap(.~Year, ncol = 1) +
#   scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
#   theme_bw() +
#   # theme(axis.text.x = element_text(angle = 90)) +
#   labs(x = "Date", y = "Percent Error [%]", color = "Model") +
#   scale_color_manual(name = "Model", 
#                      # levels = c("No Model", "2-week Validation", "Full Time Series"),
#                      values = c("No Model" = brewer.pal(8,"Dark2")[1], "2-week Validation" = brewer.pal(8,"Dark2")[6], "Full Time Series" = brewer.pal(8,"Dark2")[3]),
#                      labels = c("No Model", "2-week Validation", "Full Time Series"),
#                      guide = "legend") +
#   # scale_linetype(labels = c("2021", "2022", "2023"), guide = "legend") +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold"))
# 
# 
# 
# 
# ggplot(data = cross.val.long) +
#   geom_line(data = combined.df.overlap, aes(x = Date.Time, y = aou), lwd = 1, alpha = 0, color = "white") +
#   geom_bar(aes(x = Date.Time, y = abs(my.perc.error)*100, fill = factor(model, levels = c("2-week Validation", "Full Time Series"))), size = 2, stat = "identity") +
#   geom_hline(yintercept = 0) + 
#   # geom_line(aes(x = Date, y = my.RMSE, color = model), size = 1) +
#   # geom_point(aes(x = Date.Time, y = abs(my.perc.error)*100, color = factor(model, levels = c("2-week Validation", "Full Time Series"))), size = 2) +
#   # geom_line(aes(x = Date.Time, y = abs(my.perc.error)*100, color = factor(model, levels = c("No Model", "2-week Validation", "Full Time Series"))), size = 1) +
#   # geom_point(aes(x = Date, y = my.perc.error, color = model), size = 2) +
#   # facet_wrap(.~Year, ncol = 1) +
#   # scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
#   theme_bw() +
#   # theme(axis.text.x = element_text(angle = 90)) +
#   labs(x = "Date", y = "Percent Error [%]", color = "Model") +
#   scale_fill_manual(name = "Model", 
#                     # levels = c("No Model", "2-week Validation", "Full Time Series"),
#                     values = c("2-week Validation" = brewer.pal(8,"Dark2")[6], "Full Time Series" = brewer.pal(8,"Dark2")[3]),
#                     labels = c("2-week Validation", "Full Time Series"),
#                     guide = "legend") +
#   # scale_linetype(labels = c("2021", "2022", "2023"), guide = "legend") +
#   # theme(panel.grid = element_blank()) +
#   theme(legend.position = "none") +
#   theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold"))
# 
# 
# plot.4a <- ggplot(data = cross.val.long) +
#   geom_histogram(aes(abs(my.perc.error*100)), color = "black", fill = "#e6ab02") +
#   theme_bw() +
#   labs(x = "Percent Error [%]", y = "Frequency") +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) 
# plot.4a
# 
# combined.df$error.c <- abs(combined.df$BOU - combined.df$BOU.estimated)/abs(combined.df$BOU)
# 
# plot.4c <- ggplot(data = combined.df) +
#   geom_histogram(aes(abs(error.c)), color = "black", fill = "#e6ab02", binwidth = 5) +
#   theme_bw() +
#   labs(x = "Percent Error [%]", y = "Frequency") +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) 
# 
# # ---- PCA of predictors to compare train and test data ----
# 
# my.date.labels <- combined.df.daily$Date
# env.pred.matrix <- as.matrix(combined.df.daily[,predictors])
# env.pred.matrix <- env.pred.matrix[,-9]
# 
# PCA <- rda(env.pred.matrix, scale = FALSE)
# 
# barplot(as.vector(PCA$CA$eig)/sum(PCA$CA$eig))
# sum((as.vector(PCA$CA$eig)/sum(PCA$CA$eig))[1:2])
# smry <- summary(PCA)
# df1  <- data.frame(smry$sites[,1:2])# PC1 and PC2
# df1$Date <- my.date.labels
# df2  <- data.frame(smry$species[,1:2])     # loadings for PC1 and PC2
# df1$data.set <- "train"
# df1$data.set[which(df1$Date.Time >= date1 & df1$Date.Time <= date2)] <- "test"
# df1$data.set <- factor(df1$data.set, levels = c("train", "test"))
# 
# # df2 <- df2[which((abs(df2$PC1) + abs(df2$PC2)) >= 0.3),]
# 
# ggplot() + 
#   # geom_point(aes(color = data.set, alpha = data.set)) +
#   geom_point(data = df1[which(df1$data.set == "train"),], aes(x=PC1, y=PC2, color = "Training"), alpha = 0.3) +
#   geom_point(data = df1[which(df1$data.set == "test"),], aes(x=PC1, y=PC2, color = "Validation"), alpha = 0.8) +
#   geom_hline(yintercept=0, linetype="dotted") +
#   geom_vline(xintercept=0, linetype="dotted") +
#   coord_fixed() +
#   # geom_segment(data=df2, aes(x=0, xend=PC1, y=0, yend=PC2),
#   #              color="black", arrow=arrow(length=unit(0.01,"npc"))) +
#   # geom_text(data=df2,
#   #           aes(x=PC1,y=PC2,label=rownames(df2)),
#   #               # hjust= (1*sign(PC1)), vjust= (-0.5*sign(PC1))),
#   #           color="black", size=4) +
#   xlim(c(-max(c(df1$PC1, df1$PC2)), max(c(df1$PC1, df1$PC2)))) + ylim(c(-max(c(df1$PC1, df1$PC2)), max(c(df1$PC1, df1$PC2)))) +
#   # xlim(c(-max(c(df2$PC1, df2$PC2)), max(c(df2$PC1, df2$PC2)))) + ylim(c(-max(c(df2$PC1, df2$PC2)), max(c(df2$PC1, df2$PC2)))) +
#   theme_bw() +
#   # scale_alpha_manual(name = "Dataset", labels = factor(c("Training", "Validation"), levels = c("Training", "Validation")), guide = "legend", values = c(0.2, 0.8)) +
#   scale_color_manual(name = "Dataset", labels = factor(c("Training", "Validation"), levels = c("Training", "Validation")), guide = "legend", values = c("black", brewer.pal(n = 8, name = "Dark2")[6])) +
#   theme(panel.grid = element_blank()) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) 
# 
# 
# df1 <- merge(df1, combined.df.daily, by = "Date")
# 
# summary(lm(df1$aou~df1$PC1))
# summary(lm(df1$aou~df1$PC2))
# 
# summary(lm(df1$temperature.sccoos~df1$PC1))
# summary(lm(df1$temperature.sccoos~df1$PC2))
# 
# summary(lm(df1$salinity~df1$PC1))
# summary(lm(df1$salinity~df1$PC2))
# 
# summary(lm(df1$DPD~df1$PC1))
# summary(lm(df1$DPD~df1$PC2))
# 
# summary(lm(df1$pressure~df1$PC1))
# summary(lm(df1$pressure~df1$PC2))
# 
# summary(lm(df1$WSPD~df1$PC1))
# summary(lm(df1$WSPD~df1$PC2))
# 
# summary(lm(df1$GST~df1$PC1))
# summary(lm(df1$GST~df1$PC2))
# 
# summary(lm(df1$WVHT~df1$PC1))
# summary(lm(df1$WVHT~df1$PC2))
# 
# 
# 
# 
# df2$PC1 <- df2$PC1 / max(abs(df2$PC1)) * max(df1$PC1)
# df2$PC2 <- df2$PC2 / max(abs(df2$PC2)) * max(df1$PC2)
# df2 <- df2[which((abs(df2$PC1) + abs(df2$PC2)) >= 0.5),]
# rownames(df2) <- c("AOU", "Water Temperature", "DPD")
# 
# ggplot(df1, aes(x=PC1, y=PC2)) + 
#   geom_point(aes(color = BOU), alpha = 0.7) +
#   geom_hline(yintercept=0, linetype="dotted") +
#   geom_vline(xintercept=0, linetype="dotted") +
#   coord_fixed() +
#   geom_segment(data=df2, aes(x=0, xend=PC1, y=0, yend=PC2),
#                color="black", arrow=arrow(length=unit(0.01,"npc"))) +
#   geom_text(data=df2,
#             aes(x=PC1,y=PC2,label=rownames(df2)),
#             # hjust= (1*sign(PC1)), vjust= (-0.5*sign(PC1))),
#             color="black", size=4) +
#   # scale_color_manual(values = c("grey69", "blue")) +
#   # scale_alpha_manual(values = c(0.1, 0.8)) +
#   scale_color_viridis_c() +
#   # xlim(c(-max(c(df2$PC1, df2$PC2)), max(c(df2$PC1, df2$PC2)))) + ylim(c(-max(c(df2$PC1, df2$PC2)), max(c(df2$PC1, df2$PC2)))) +
#   theme_bw() +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   theme(panel.grid = element_blank())
# 
# plot.4d <- ggplot(df1, aes(x=PC1, y=PC2)) + 
#   geom_point(aes(color = BOU.estimated), alpha = 0.7) +
#   geom_hline(yintercept=0, linetype="dotted") +
#   geom_vline(xintercept=0, linetype="dotted") +
#   # coord_fixed() +
#   geom_segment(data=df2, aes(x=0, xend=PC1, y=0, yend=PC2),
#                color="black", arrow=arrow(length=unit(0.01,"npc"))) +
#   geom_text(data=df2,
#             aes(x=PC1,y=PC2,label=rownames(df2)),
#             # hjust= (1*sign(PC1)), vjust= (-0.5*sign(PC1))),
#             color="black", size=4) +
#   # scale_color_manual(values = c("grey69", "blue")) +
#   # scale_alpha_manual(values = c(0.1, 0.8)) +
#   scale_color_viridis_c() +
#   labs(color = "Estimated \n BOU") +
#   # xlim(c(-max(c(df2$PC1, df2$PC2)), max(c(df2$PC1, df2$PC2)))) + ylim(c(-max(c(df2$PC1, df2$PC2)), max(c(df2$PC1, df2$PC2)))) +
#   theme_bw() +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   theme(panel.grid = element_blank())
# 
# df1$error.c <- 100*abs(df1$BOU - df1$BOU.estimated)/abs(df1$BOU)
# df1$na.alpha <- 0.7
# df1$na.alpha[which(is.na(df1$error.c))] <- 0.3
# 
# 
# plot.4e <- ggplot(df1, aes(x=PC1, y=PC2)) + 
#   geom_point(aes(color = log(abs(error.c)), alpha = na.alpha)) +
#   geom_hline(yintercept=0, linetype="dotted") +
#   geom_vline(xintercept=0, linetype="dotted") +
#   # coord_fixed() +
#   geom_segment(data=df2, aes(x=0, xend=PC1, y=0, yend=PC2),
#                color="black", arrow=arrow(length=unit(0.01,"npc"))) +
#   geom_text(data=df2,
#             aes(x=PC1,y=PC2,label=rownames(df2)),
#             # hjust= (01*sign(PC1)), vjust= (-0.5*sign(PC1))),
#             color="black", size=4) +
#   # scale_color_manual(values = c("grey69", "blue")) +
#   # scale_alpha_manual(values = c(0.1, 0.8)) +
#   scale_color_viridis_c() +
#   # xlim(c(-max(c(df2$PC1)), max(c(df2$PC1)))) + ylim(c(-max(c(df2$PC2)), max(c(df2$PC2)))) +
#   labs(color = paste0(" Log \n Relative \n Error [%]")) +
#   theme_bw() +
#   scale_alpha(guide = "none") +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   theme(panel.grid = element_blank())
# 
# cross.val.long.to.combine <- cross.val.long[,c(1,3)]
# colnames(cross.val.long.to.combine)[1] <- "Date"
# 
# df1 <- merge(df1, cross.val.long.to.combine, by = "Date", all = T)
# 
# df1$na.alpha <- 0.7
# df1$na.alpha[which(is.na(df1$my.perc.error))] <- 0.3
# 
# plot.4b <- ggplot(df1, aes(x=PC1, y=PC2)) + 
#   geom_point(aes(color = log(abs(my.perc.error*100)), alpha = na.alpha), size = 2) +
#   geom_hline(yintercept=0, linetype="dotted") +
#   geom_vline(xintercept=0, linetype="dotted") +
#   # coord_fixed() +
#   geom_segment(data=df2, aes(x=0, xend=PC1, y=0, yend=PC2),
#                color="black", arrow=arrow(length=unit(0.01,"npc"))) +
#   geom_text(data=df2,
#             aes(x=PC1,y=PC2,label=rownames(df2)),
#             # hjust= (1*sign(PC1)), vjust= (-0.5*sign(PC1))),
#             color="black", size=4) +
#   # scale_color_manual(values = c("grey69", "blue")) +
#   # scale_alpha_manual(values = c(0.1, 0.8)) +
#   scale_color_viridis_c() +
#   # xlim(c(-max(c(df2$PC1)), max(c(df2$PC1)))) + ylim(c(-max(c(df2$PC2)), max(c(df2$PC2)))) +
#   labs(color = paste0(" Log \n Relative \n Error [%]")) +
#   theme_bw() +
#   scale_alpha(guide = "none") +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   theme(panel.grid = element_blank())
# 
# 
# # (plot.4a + plot.4c) / (plot.4b + plot.4d + plot.4e) + plot_layout(heights = c(2,4))
# 
# png("../Figures/manuscript_figures/2024-08-09_fig4.png", width = 1300, height = 500)
# 
# (plot.4a + plot.4b) 
# 
# dev.off()
# 
# png("../Figures/manuscript_figures/2024-08-09_fig6.png", width = 1300, height = 500)
# 
# plot.4d + (plot.4c / plot.4e) 
# 
# dev.off()
# 
# 
# summary(lm(df1$my.perc.error ~ df1$temperature.sccoos))
# summary(lm(df1$my.perc.error ~ df1$salinity))
# summary(lm(df1$my.perc.error ~ df1$pressure))
# summary(lm(df1$my.perc.error ~ df1$DPD))
# summary(lm(df1$my.perc.error ~ df1$WSPD))
# summary(lm(df1$my.perc.error ~ df1$GST))
# summary(lm(df1$my.perc.error ~ df1$WVHT))
# summary(lm(df1$my.perc.error ~ df1$BOU))
# 
# 
# 
# 
# 
# # ---- MISC ----
# 
# try.it <- cross.val.long
# colnames(try.it)[1] <- "Date"
# try.it <- try.it[,-4]
# try.it <- merge(x = try.it, y = combined.df.daily, by = "Date", all.x = T, all.y = F)
# 
# summary(lm(try.it$my.perc.error*100~try.it$BOU))
# plot(try.it$BOU, try.it$my.perc.error*100)
# abline(a = summary(lm(try.it$my.perc.error*100~try.it$BOU))$coefficients[1,1], summary(lm(try.it$my.perc.error*100~try.it$BOU))$coefficients[2,1])
# 
# 
# try.it$signs.match <- "TRUE"
# try.it$signs.match[which(sign(try.it$BOU) != sign(try.it$BOU.estimated))] <- "FALSE"
# 
# ggplot(data = try.it) +
#   # geom_line(aes(x = Date, y = my.RMSE, color = model), size = 1) +
#   geom_boxplot(aes(color = factor(signs.match, levels = c("TRUE", "FALSE")), 
#                    y = abs(my.perc.error)*100, 
#                    x = factor(model, levels = c("2-week Validation", "Full Time Series")), 
#   )) +
#   # geom_point(aes(x = Date, y = my.perc.error, color = model), size = 2) +
#   # facet_wrap(.~Year, ncol = 1) +
#   # scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
#   theme_bw() +
#   # theme(axis.text.x = element_text(angle = 90)) +
#   labs(x = "Model", y = "Percent Error [%]", color = "Same Sign?") +
#   scale_x_discrete(labels = c("2-week Validation", "Full Time Series")) +
#   # scale_color_manual(values = c("2-week Validation" = brewer.pal(8,"Dark2")[6], "Full Time Series" = brewer.pal(8,"Dark2")[3]),
#   #                    labels = c("2-week Validation", "Full Time Series"),
#   #                    guide = "legend") +
#   scale_color_manual(values = c("TRUE" = brewer.pal(8,"Dark2")[6], "FALSE" = brewer.pal(8,"Dark2")[3]),
#                      labels = c("TRUE", "FALSE"),
#                      guide = "legend") +
#   # scale_linetype(labels = c("2021", "2022", "2023"), guide = "legend") +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold"))
# 
# my.aov <- aov(try.it$my.perc.error~try.it$signs.match)
# summary(my.aov)
# 
# 
# 
# cross.val.long$Season <- "Winter"
# cross.val.long$Season[which(month(cross.val.long$Date.Time) >= 4 & month(cross.val.long$Date.Time) <= 9)] <- "Summer"
# 
# 
# ggplot(data = cross.val.long) +
#   # geom_line(aes(x = Date, y = my.RMSE, color = model), size = 1) +
#   geom_boxplot(aes(color = factor(Season, levels = c("Summer", "Winter")), 
#                    y = abs(my.perc.error)*100, 
#                    x = factor(model, levels = c("2-week Validation", "Full Time Series")), 
#   )) +
#   # geom_point(aes(x = Date, y = my.perc.error, color = model), size = 2) +
#   # facet_wrap(.~Year, ncol = 1) +
#   # scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
#   theme_bw() +
#   # theme(axis.text.x = element_text(angle = 90)) +
#   labs(x = "Model", y = "Percent Error [%]", color = "Season") +
#   scale_x_discrete(labels = c("2-week Validation", "Full Time Series")) +
#   # scale_color_manual(values = c("2-week Validation" = brewer.pal(8,"Dark2")[6], "Full Time Series" = brewer.pal(8,"Dark2")[3]),
#   #                    labels = c("2-week Validation", "Full Time Series"),
#   #                    guide = "legend") +
#   scale_color_manual(values = c("Summer" = brewer.pal(8,"Dark2")[6], "Winter" = brewer.pal(8,"Dark2")[3]),
#                      labels = c("Summer", "Winter"),
#                      guide = "legend") +
#   # scale_linetype(labels = c("2021", "2022", "2023"), guide = "legend") +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold"))
# 
# combined.df.select$Season <- "Winter"
# combined.df.select$Season[which(month(combined.df.select$Date.Time) >= 4 & month(combined.df.select$Date.Time) <= 9 )] <- "Summer"
# nrow(combined.df.select[which(combined.df.select$Season == "Summer"),])
# nrow(combined.df.select[which(combined.df.select$Season == "Winter"),])
# nrow(combined.df.select[which(combined.df.select$Season == "Winter"),])/nrow(combined.df.select[which(combined.df.select$Season == "Summer"),])
# 
# two.week.df <- cross.val.long[which(cross.val.long$model == "2-week Validation"),]
# # full.ts.df <- cross.val.long[which(cross.val.long$model == "Full Time Series"),]
# 
# t.test(x = two.week.df$my.perc.error[which(two.week.df$Season == "Summer")], y = two.week.df$my.perc.error[which(two.week.df$Season == "Winter")])
# # t.test(x = full.ts.df$my.perc.error[which(full.ts.df$Season == "Summer")], y = full.ts.df$my.perc.error[which(full.ts.df$Season == "Winter")])
# 
# 
# 
# 
# 
# 
# my.aov <- aov(cross.val.long$my.perc.error~cross.val.long$model)
# TukeyHSD(my.aov)
# 
# 
# 
# df1$Date.Time <- parse_date_time(paste(year(df1$Date.Time), month(df1$Date.Time), day(df1$Date.Time), sep = "-"), orders = "Ymd")
# df1 <- df1 %>% group_by(Date.Time) %>% summarize(PC1 = mean(PC1), PC2 = mean(PC2))
# 
# try.it <- merge(df1, cross.val.long, by = "Date.Time")
# 
# ggplot(data = try.it[which(try.it$model == "RMSE_model_hyper_full"),], aes(x = PC1, y = PC2)) +
#   geom_point(aes(color = my.perc.error), size = 2) +
#   geom_hline(yintercept=0, linetype="dotted") +
#   geom_vline(xintercept=0, linetype="dotted") +
#   coord_fixed() +
#   geom_segment(data=df2, aes(x=0, xend=PC1/10, y=0, yend=PC2/10),
#                color="black", arrow=arrow(length=unit(0.01,"npc"))) +
#   geom_text(data=df2,
#             aes(x=PC1/10,y=PC2/10,label=rownames(df2),
#                 hjust= (1*sign(PC1)), vjust= (-0.5*sign(PC1))),
#             color="black", size=4) +
#   # scale_color_manual(values = c("grey69", "blue")) +
#   scale_color_viridis_c() +
#   # scale_alpha_manual(values = c(0.1, 0.8)) +
#   xlim(c(-max(c(df1$PC1, df1$PC2)), max(c(df1$PC1, df1$PC2)))) + ylim(c(-max(c(df1$PC1, df1$PC2)), max(c(df1$PC1, df1$PC2)))) +
#   theme_bw()
# 
# 
# combined.df.select2 <- combined.df.select
# combined.df.select2$Date.Time <- parse_date_time(paste(year(combined.df.select$Date.Time), month(combined.df.select$Date.Time), day(combined.df.select$Date.Time), sep = "-"), orders = "Ymd")
# combined.df.select2 <- combined.df.select2 %>% group_by(Date.Time) %>% summarize_all(mean)
# 
# try.it2 <- merge(try.it, combined.df.select2, by = "Date.Time", all.x = T, ally. = F)
# try.it2 <- try.it2 %>% pivot_wider(names_from = model, values_from = my.perc.error)
# try.it2$model.diff.hyper <- try.it2$avg_perc_error_no_model - try.it2$avg_perc_error_hyper  # higher is better
# try.it2$model.diff.hyper.full <- try.it2$avg_perc_error_no_model - try.it2$avg_perc_error_hyper_full  # higher is better
# 
# 
# summary(lm(formula = try.it2$model.diff.hyper~try.it2$BOU))
# summary(lm(formula = try.it2$model.diff.hyper.full~try.it2$BOU))
# 
# plot.3c <- ggplot(data = try.it2) +
#   geom_point(aes(x = BOU, y = model.diff.hyper, color = "2-week Validation Period - No Model"), size = 2) +
#   geom_smooth(aes(x = BOU, y = model.diff.hyper, color = "2-week Validation Period - No Model"), se = F, method = "lm", linewidth = 2) +
#   geom_point(aes(x = BOU, y = model.diff.hyper.full, color = "Full Time Series Model - No Model"), size = 2) +
#   geom_smooth(aes(x = BOU, y = model.diff.hyper.full, color = "Full Time Series Model - No Model"), se = F, method = "lm", linewidth = 2) +
#   theme_bw() +
#   labs(x = "[O2]bio [μM]", y = "RMSE Difference [μM]", color = "") +
#   scale_color_manual(values = c(col5.other, "black"), labels = c("2-week Validation Period - No Model", "Full Time Series Model - No Model"), guide = "legend") +
#   # scale_linetype(labels = c("2021", "2022", "2023"), guide = "legend") +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   theme(axis.title.y = element_blank(), axis.text.y = element_blank())
# 
# summary(lm(formula = try.it2$model.diff.hyper~try.it2$aou))
# summary(lm(formula = try.it2$model.diff.hyper.full~try.it2$aou))
# 
# plot.3b <- ggplot(data = try.it2) +
#   geom_point(aes(x = aou, y = model.diff.hyper, color = "2-week Validation Period - No Model"), size = 2) +
#   geom_smooth(aes(x = aou, y = model.diff.hyper, color = "2-week Validation Period - No Model"), se = F, method = "lm", linewidth = 2) +
#   geom_point(aes(x = aou, y = model.diff.hyper.full, color = "Full Time Series Model - No Model"), size = 2) +
#   geom_smooth(aes(x = aou, y = model.diff.hyper.full, color = "Full Time Series Model - No Model"), se = F, method = "lm", linewidth = 2) +
#   theme_bw() +
#   labs(x = "aou", y = "RMSE Difference [μM]", color = "") +
#   scale_color_manual(values = c(col5.other, "black"), labels = c("2-week Validation Period - No Model", "Full Time Series Model - No Model"), guide = "legend") +
#   # scale_linetype(labels = c("2021", "2022", "2023"), guide = "legend") +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   theme(legend.position = "none")
# 
# plot.3a / (plot.3b + plot.3c)
# 
# 
# 
# ggplot(data = try.it2) +
#   geom_point(aes(x = temperature.sccoos, y = BOU), color = "black", size = 2) +
#   geom_smooth(aes(x = temperature.sccoos, y = BOU), color = "black", se = F, method = "lm", linewidth = 2) +
#   theme_bw() +
#   labs(x = "Water Temperature [C]", y = "[O2]bio [μM]") +
#   # scale_linetype(labels = c("2021", "2022", "2023"), guide = "legend") +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) 
# 
# summary(lm(formula = try.it2$temperature.sccoos~try.it2$BOU))
# 
# ggplot(data = combined.df.select) +
#   geom_point(aes(x = temperature.sccoos, y = BOU), color = "black", size = 2) +
#   geom_smooth(aes(x = temperature.sccoos, y = BOU), color = "red", se = F, method = "lm", linewidth = 2) +
#   theme_bw() +
#   labs(x = "Water Temperature [C]", y = "[O2]bio [μM]") +
#   # scale_linetype(labels = c("2021", "2022", "2023"), guide = "legend") +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) 
# 
# summary(lm(formula = combined.df.select$temperature.sccoos~combined.df.select$BOU))
# 
# 
# hist(try.it2$temperature.sccoos)
# hist(try.it2$BOU)
# hist(try.it2$RMSE_no_model)
# hist(try.it2$RMSE_model_hyper)
# hist(try.it2$RMSE_model_hyper_full)
# 
# 
# 
# 
# 
# ggplot(data = try.it2) +
#   geom_point(aes(x = salinity, y = model.diff.hyper, color = "2-week Validation Period - No Model"), size = 2) +
#   geom_smooth(aes(x = salinity, y = model.diff.hyper, color = "2-week Validation Period - No Model"), se = F, method = "lm", linewidth = 2) +
#   geom_point(aes(x = salinity, y = model.diff.hyper.full, color = "Full Time Series Model - No Model"), size = 2) +
#   geom_smooth(aes(x = salinity, y = model.diff.hyper.full, color = "Full Time Series Model - No Model"), se = F, method = "lm", linewidth = 2) +
#   theme_bw() +
#   labs(x = "Salinity [psu]", y = "RMSE Difference [μM]", color = "") +
#   scale_color_manual(values = c(col5.other, "black"), labels = c("2-week Validation Period - No Model", "Full Time Series Model - No Model"), guide = "legend") +
#   # scale_linetype(labels = c("2021", "2022", "2023"), guide = "legend") +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) 
# 
# 
# 
# # ---- smoothing ----
# 
# date.start <- full.predictors$Date.Time[1] - weeks(1)
# date.end <- full.predictors$Date.Time[nrow(full.predictors)] + weeks(1)
# test <- seq(from = date.start, to = date.end, by = "hours")
# 
# full.predictors.with.BOU.est <- full.predictors
# full.predictors.with.BOU.est$aou_corrected <- full.BOU.estrected
# 
# my.rolling.avgs.df <- full.predictors.with.BOU.est
# 
# # clunky but it works
# for(r in 1:nrow(full.predictors.with.BOU.est)){
#   
#   temp.df <- full.predictors.with.BOU.est[which(full.predictors.with.BOU.est$Date.Time >= (full.predictors.with.BOU.est$Date.Time[r] - weeks(1)) & full.predictors.with.BOU.est$Date.Time <= (full.predictors.with.BOU.est$Date.Time[r] + weeks(1))),]
#   
#   my.rolling.avgs.df[r,-9] <- colMeans(temp.df[-9])
#   
# }
# 
# saveRDS(my.rolling.avgs.df)
# 
# 
# sccoos.dates <- readRDS("../16S_sccoos/2023-04-28_sccoos_dates.rds")
# 
# ggplot() +
#   geom_line(aes(x = full.predictors$Date.Time, y = full.BOU.estrected), color = col3.BOU_est, lwd = 1, alpha = 0.7) +
#   geom_line(data = my.rolling.avgs.df, aes(x = Date.Time, y = aou_corrected), color = "black", linewidth = 1) +
#   # geom_line(data = my.rolling.avgs.df, aes(x = Date.Time, y = aou), color = col1) +
#   # geom_line(data = my.rolling.avgs.df, aes(x = Date.Time, y = BOU), color = col3) +
#   geom_hline(aes(yintercept = 0), alpha = 0.2) +
#   labs(x = "Date", y = expression("Corrected aou  ["*mu*"M]")) +
#   # ylim(c(-100,100)) +
#   theme_bw() +
#   scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
#   geom_jitter(aes(x = parse_date_time(sccoos.dates, orders = "Ymd"), y = 140), height = 10, alpha = 0.7) +
#   theme(panel.grid = element_blank()) +
#   theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12)) 
# 
# 
# 
# # ---- cartoon model ----
# 
# ## This is a toy decision tree that is useful for talking to an audience about
# ## gradient boosted trees, etc.
# 
# dummy.predictors <- c("BOU",
#                       "temperature.sccoos.zoo",
#                       'salinity',
#                       "pressure",
#                       "WSPD",
#                       "GST",
#                       "WVHT",
#                       "DPD")
# library(tree)
# 
# temp <- tree(BOU ~ .,
#              data = na.omit(combined.df[,predictors]))
# plot(temp)
# text(temp)
# 
# 
# #### export data ####
# 
# saveRDS(full.aou.df, "full_aou_df.rds")
# saveRDS(full.predictors, "full_predictors_df.rds")
# saveRDS(mims.hourly, "mims_hourly.rds")
# 
