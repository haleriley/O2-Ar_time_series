temp <- tempfile()
download.file("https://www.polarmicrobes.org/ecoobs/CTD_data_vol1.csv.gz", temp)
my.file <- read.csv(gzfile(temp), as.is = TRUE)
unlink(temp)

# temp <- tempfile()
# download.file("https://www.polarmicrobes.org/CTD_data_vol2.csv.gz", temp)
# my.file2 <- read.csv(gzfile(temp), as.is = TRUE)
# unlink(temp)

library(tidyverse)
library(lubridate)


my.file$X <- parse_date_time(substr(my.file$X, 1,19), orders = "Ymd HMS")

# download MIMS2.1 data
# temp <- tempfile()
# download.file("https://polarmicrobes.org/ecoobs/o2bio_vol2.1.csv", temp) # loads in current MIMS vol2 data
# my.file <- read.csv(gzfile(temp), as.is = TRUE)
# unlink(temp)
# 
# mims2.1 <- my.file
# # mims2 <- my.file[,c("date_time", "O2", "Temperature..ITS.90.deg.C.")]
# mims2.1$date_time <- parse_date_time(paste(year(mims2.1$date_time), "-", month(mims2.1$date_time), "-", day(mims2.1$date_time), " ", hour(mims2.1$date_time), sep = ""), orders = "Ymd H")
# # mims2 <- my.file[,c("date_time", "N2", "O2", "Ar", "Inlet.Temperature")]
# # mims2$date_time <- parse_date_time(paste(year(mims2$time), "-", month(mims2$time), "-", day(mims2$time), " ", hour(mims2$time), sep = ""), orders = "Ymd H")
# # mims2 <- mims2[,-1]
# mims2.1 <- mims2.1 %>% group_by(date_time) %>% summarize_all(mean)
# 
# my.file <- mims2.1
# colnames(my.file)[1] <- "X"

# ---- FULL TIME SERIES ----

# temp full
ggplot(my.file) +
  geom_line(aes(x = X, y = Temperature..ITS.90.deg.C.)) +
  theme_bw() +
  scale_x_datetime(date_breaks = "2 months", date_labels = "%b %Y")

# salinity full
ggplot(my.file) +
  geom_line(aes(x = X, y = Salinity..PSU.)) +
  theme_bw() +
  scale_x_datetime(date_breaks = "2 months", date_labels = "%b %Y")

# cond full
ggplot(my.file) +
  geom_line(aes(x = X, y = Conductivity..mS.cm.)) +
  theme_bw() +
  scale_x_datetime(date_breaks = "2 months", date_labels = "%b %Y")

# O2 full 
ggplot(my.file) +
  geom_line(aes(x = X, y = Oxygen..umol.l.)) +
  theme_bw() +
  scale_x_datetime(date_breaks = "2 months", date_labels = "%b %Y")

# fluor full
ggplot(my.file) +
  geom_line(aes(x = X, y = Fluorescence..mg.m.3.)) +
  theme_bw() +
  scale_x_datetime(date_breaks = "2 months", date_labels = "%b %Y")


# ---- month two weeks ----

last.date <- my.file$X[length(my.file$X)]
one.month.ago <- last.date - (2*30*24*60*60)
oma.index <- which(my.file$X >= one.month.ago)[1]

my.file.past.month <- my.file[c(oma.index:nrow(my.file)),]

# temp month
ggplot(my.file.past.month) +
  geom_line(aes(x = X, y = Temperature..ITS.90.deg.C.)) +
  theme_bw() +
  scale_x_datetime(date_breaks = "3 days", date_labels = "%d %b")

# salinity month
ggplot(my.file.past.month) +
  geom_line(aes(x = X, y = Salinity..PSU.)) +
  theme_bw() +
  scale_x_datetime(date_breaks = "3 days", date_labels = "%d %b")

# cond month
ggplot(my.file.past.month) +
  geom_line(aes(x = X, y = Conductivity..mS.cm.)) +
  theme_bw() +
  scale_x_datetime(date_breaks = "3 days", date_labels = "%d %b")

# O2 month
ggplot(my.file.past.month) +
  geom_line(aes(x = X, y = Oxygen..umol.l.)) +
  theme_bw() +
  scale_x_datetime(date_breaks = "3 days", date_labels = "%d %b")

# fluor month
ggplot(my.file.past.month) +
  geom_line(aes(x = X, y = Fluorescence..mg.m.3.)) +
  theme_bw() +
  scale_x_datetime(date_breaks = "3 days", date_labels = "%d %b")

