# 2023-08-21
# AOP Model Maximum Error 
# RJH

getwd()
setwd('C://Users/haler/Documents/PhD-Bowman/O2-Ar_time_series/R_Data/')
set.seed(1234)


# ---- library ----

library(tidyverse)
library(lubridate)
library(rLakeAnalyzer)

# ---- read in data ----

# combined.df <- readRDS(file = "2023-08-08_combined_env_data_hourly.rds")
final.df <- readRDS("2023-08-08_aop_cor_df.rds")

# ---- combine to make final dataset ----

# combo <- merge(final.df, combined.df, by = "Date.Time", all.x = T, all.y = F)
# combo <- combo[which(is.na(combo$o2_bio) == F),]
combo <- final.df

# ---- calculate final error ----

# just to make sure matches other calculations
# sqrt(mean((combo$o2_bio - combo$aop.corrected)^2)) # looks right

combo$delta.final <- combo$aop - combo$aop.corrected

summary(abs(combo$delta.final))


abs(combo$delta.final[head(order(abs(combo$delta.final), decreasing = T), n = round(0.01*nrow(combo)))])

combo$max.error.Y.N <- "N"
combo$max.error.Y.N[which(abs(combo$delta.final) > 100)] <- "Y"


ggplot(data = combo) +
  # geom_line(aes(x = Date.Time, y = abs(combo$delta.final))) +
  geom_point(aes(x = Date.Time, y = abs(combo$delta.final), color = max.error.Y.N)) +
  theme_bw()


high.error <- combo[which(combo$max.error.Y.N == "Y"),]















