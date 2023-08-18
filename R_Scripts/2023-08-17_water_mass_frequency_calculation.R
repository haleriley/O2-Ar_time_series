# water mass frequency calibration #
# 2023-08-17
# RJH

getwd()
setwd('C://Users/haler/Documents/PhD-Bowman/O2-Ar_time_series/R_Data/')
set.seed(1234)


# ---- library ----

library(tidyverse)
library(lubridate)
library(rLakeAnalyzer)

# ---- read in data ----

# final.df <- readRDS(file = "2023-08-08_combined_env_data_daily.rds")
final.df <- readRDS("2023-08-08_aop_cor_df.rds")


final.df$Date.Time <- parse_date_time(paste(year(final.df$Date.Time), month(final.df$Date.Time), day(final.df$Date.Time), sep = "-"), orders = "Ymd")
final.df <- final.df %>% group_by(Date.Time) %>% summarize_all(mean)

# ---- calculate density ----

final.df$w.dens <- water.density(wtr = final.df$temperature.sccoos, sal = final.df$salinity)


# ---- view in T-S space ----

ggplot(data = final.df) +
  geom_point(aes(x = temperature.sccoos, y = salinity, color = month(Date.Time))) +
  facet_wrap(.~year(Date.Time)) +
  theme_bw()


# ---- calculate change in temp and salinity ----

final.df$temp.diff <- c(NA, diff(final.df$temperature.sccoos)/as.numeric(diff(final.df$Date.Time)))
final.df$sal.diff <- c(NA, diff(final.df$salinity)/as.numeric(diff(final.df$Date.Time)))

summary(final.df$temp.diff)
hist(final.df$temp.diff)
sd(na.omit(final.df$temp.diff))
sd.temp <- sd(na.omit(final.df$temp.diff))


summary(final.df$sal.diff)
hist(final.df$sal.diff)
sd(na.omit(final.df$sal.diff))
sd.sal <- sd(na.omit(final.df$sal.diff))


test1 <- final.df[which(abs(final.df$temp.diff) > sd.temp & abs(final.df$sal.diff) > sd.sal),]

Date.Time.diff1 <- diff(test1$Date.Time)
mean(Date.Time.diff1)
median(Date.Time.diff1)


# ---- aop change frequency ----

final.df$aop.cor.diff <- c(NA, diff(final.df$aop.corrected)/as.numeric(diff(final.df$Date.Time)))

summary(final.df$aop.cor.diff)
hist(final.df$aop.cor.diff)
sd(na.omit(final.df$aop.cor.diff))
sd.aop <- sd(na.omit(final.df$aop.cor.diff))

test.aop <- final.df[which(abs(final.df$aop.cor.diff) > sd.aop),]

Date.Time.diff.aop <- diff(test.aop$Date.Time)
mean(Date.Time.diff.aop)
median(Date.Time.diff.aop)


# ---- compare ----

summary(lm(abs(final.df$temp.diff)~abs(final.df$sal.diff)))
plot(abs(final.df$temp.diff), abs(final.df$sal.diff))

summary(lm(abs(final.df$temp.diff)~abs(final.df$aop.cor.diff)))
plot(abs(final.df$temp.diff), abs(final.df$aop.cor.diff))

summary(lm(abs(final.df$sal.diff)~abs(final.df$aop.cor.diff)))
plot(abs(final.df$sal.diff), abs(final.df$aop.cor.diff))












