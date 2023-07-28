# 2022-10-21
# AOU temporal patterns
# RJH

getwd()
setwd('C://Users/haler/Documents/PhD-Bowman/MIMS-miniDOT_O2-Ar_Study/R_Data/')

aou.df <- read.csv(file = "2023-04-28_corrected_aou.csv")


col1 <- "#2274A5"
col2 <- "#F75C03"
col3 <- "#F1C40F"
col4 <- "#D90368"
col5 <- "#00CC66"

# ---- library ----

library(tidyverse)
library(lubridate)
library(zoo)
library(plotly)
library(patchwork)
library(mgcv)
library(nls2)

# ---- format data ----

aou.df$dates <- parse_date_time(aou.df$Date.Time, orders = "Ymd HMS")

# aou.df$aou.cor.rolling.avg <- c(rollmean(x = aou.df$aou.corrected, k = 48), rep(NA, times = 47))

aou.df$NCP.status <- "Net Autotrophic"
aou.df$NCP.status[which(aou.df$aou.corrected < 0)] <- "Net Heterotrophic"


aou.df$delta <- aou.df$aou - aou.df$aou.corrected

# aou.df$model.performance <- aou.df$o2_bio - aou.df$aou.corrected
# aou.df$no.model.performance <- aou.df$o2_bio - aou.df$aou


# ---- calculate rolling average ----

aou.df <- aou.df[,-which(substr(colnames(aou.df), start = 1, stop = 5) == "Group")]

try.it <- data.frame(matrix(NA, nrow = nrow(aou.df), ncol = 26))
colnames(try.it) <- c("dates", colnames(aou.df)[c(4:length(aou.df)-1)])

for(t in 1:nrow(aou.df)){
  
  my.date <- aou.df$dates[t]
  
  my.lower.date <- aou.df$dates[t] - 12*60*60
  my.upper.date <- aou.df$dates[t] +12*60*60
  
  my.df <- aou.df[which(aou.df$dates >= my.lower.date & aou.df$dates <= my.upper.date),]
  
  test <- colMeans(my.df[,c(4:length(my.df)-1)])
  
  try.it[t,1] <- as.character(my.date)
  
  try.it[t,2:26] <- test
  
}

aou.df <- try.it

aou.df$NCP.status <- "Net Autotrophic"
aou.df$NCP.status[which(aou.df$aou.corrected < 0)] <- "Net Heterotrophic"

aou.df$dates <- parse_date_time(aou.df$dates, orders = "Ymd HMS")

# ---- plot ----

ggplot(data = aou.df) +
  geom_line(aes(x = dates, y = aou, color = NCP.status)) +
  scale_color_manual(values = c("blue", "red")) +
  geom_smooth(aes(x = dates, y = aou), color = "black", size = 2) +
  geom_hline(yintercept = 0, color = "black") +
  labs(x = "Date", y = "AOU", color = "NCP Status") +
  theme_bw()

ggplot(data = aou.df) +
  geom_line(aes(x = dates, y = aou.corrected)) +
  scale_color_manual(values = c(col5, col2)) +
  # geom_smooth(aes(x = dates, y = aou.corrected), color = "black", size = 2) + #default chooses method = gam() --> generalized additive model
  # geom_smooth(aes(x = dates, y = aou.corrected), color = "black", size = 1, method = lm, linetype = "dashed") +
  geom_hline(yintercept = 0, color = "black") +
  # geom_hline(yintercept = mean(aou.df$aou.corrected), color = "black", size = 2) +
  labs(x = "Date", y = "Corrected AOP", color = "NCP Status") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12)) 

mean(aou.df$aou.corrected)

ggplot(data = aou.df) +
  geom_line(aes(x = dates, y = aou.corrected, color = NCP.status)) +
  scale_color_manual(values = c(col5, col2)) +
  geom_smooth(aes(x = dates, y = aou.corrected), color = "black", size = 2, method = "gam") + #default chooses method = gam() --> generalized additive model
  # geom_smooth(aes(x = dates, y = aou.corrected), color = "black", size = 1, method = lm, linetype = "dashed") +
  geom_hline(yintercept = 0, color = "black") +
  labs(x = "Date", y = "Corrected AOU", color = "NCP Status") +
  theme_bw()

ggplot(data = aou.df) +
  geom_line(aes(x = dates, y = aou.corrected), color = col5) +
  # scale_color_manual(values = c(col5, col2)) +
  geom_smooth(aes(x = dates, y = aou.corrected), color = "black", size = 2, method = "gam") + #default chooses method = gam() --> generalized additive model
  # geom_smooth(aes(x = dates, y = aou.corrected), color = "black", size = 1, method = lm, linetype = "dashed") +
  geom_hline(yintercept = 0, color = "black") +
  labs(x = "Date", y = expression("Corrected AOP  [ "*mu*"M]"), color = "NCP Status") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12)) 

sccoos.dates <- readRDS("../16S_sccoos/2023-04-28_sccoos_dates.rds")
sccoos.dates <- read_table("../16S_sccoos/R_Data/2023-04-28_SCCOOS_dates.txt", col_names = F)
test <- paste(sccoos.dates$X1, "000000")

ggplot() +
  geom_line(data = aou.df, aes(x = dates, y = aou.corrected), color = "grey69") +
  # scale_color_manual(values = c(col5, col2)) +
  # geom_smooth(aes(x = dates, y = aou.corrected), color = "black", size = 2, method = "gam") + #default chooses method = gam() --> generalized additive model
  # geom_smooth(aes(x = dates, y = aou.corrected), color = "black", size = 1, method = lm, linetype = "dashed") +
  geom_jitter(aes(x = parse_date_time(test, orders = "ymd HMS"), y = 150), color = "black", width = 0, height = 20) +
  # geom_hline(yintercept = 0, color = "black") +
  labs(x = "Date", y = expression("Corrected AOP  [ "*mu*"M]"), color = "NCP Status") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12)) 





summary(lm(formula = aou.corrected~dates, data = aou.df))
summary(gam(formula = aou.corrected~dates, data = aou.df))


# temperature
ggplot(data = aou.df) +
  geom_line(aes(x = dates, y = temperature.sccoos.zoo), color = "red") +
  # scale_color_manual(values = c("blue", "red")) +
  geom_smooth(aes(x = dates, y = temperature.sccoos.zoo), color = "black", size = 2) + #default chooses method = gam() --> generalized additive model
  geom_smooth(aes(x = dates, y = temperature.sccoos.zoo), color = "black", size = 1, method = lm, linetype = "dashed") +
  # geom_hline(yintercept = 0, color = "black") +
  labs(x = "Date", y = "Temperature (SCCOOS) [C]", color = "NCP Status") +
  theme_bw()

summary(lm(formula = temperature.sccoos.zoo~dates, data = aou.df))
summary(gam(formula = temperature.sccoos.zoo~dates, data = aou.df))


ggplot(data = aou.df) +
  geom_point(aes(x = temperature.sccoos.zoo, y = o2_bio, color = NCP.status))+
  geom_hline(yintercept = 0, color = "black") +
  labs(x = "Temperature (SCCOOS) [C]", y = "AOU (corrected)", color = "NCP Status") +
  scale_color_manual(values = c("blue", "red")) +
  geom_smooth(aes(x = temperature.sccoos.zoo, y = aou.corrected), color = "black", size = 1, method = lm) +
  theme_bw()

summary(lm(formula = aou.corrected~temperature.sccoos.zoo, data = aou.df))


ggplot(data = aou.df) +
  geom_line(aes(x = dates, y = aou.corrected)) +
  # scale_color_manual(values = c("blue", "red")) +
  geom_smooth(aes(x = dates, y = aou.corrected), color = "blue") +
  geom_hline(yintercept = 0, color = "red") +
  theme_bw()

ggplot(data = aou.df) +
  geom_line(aes(x = dates, y = aou.cor.rolling.avg)) +
  scale_color_manual(values = c("blue", "red")) +
  geom_smooth(aes(x = dates, y = aou.corrected), color = "blue") +
  geom_hline(yintercept = 0, color = "red") +
  theme_bw()


ggplot(data = aou.df) +
  geom_line(aes(x = dates, y = delta)) +
  # scale_color_manual(values = c("blue", "red")) +
  # geom_smooth(aes(x = dates, y = delta), color = "black", size = 2) +
  geom_hline(yintercept = 0, color = "black") +
  labs(x = "Date", y = "AOU - AOOU (corrected)", color = "NCP Status") +
  theme_bw()

ggplot(data = aou.df) +
  geom_line(aes(x = dates, y = abs(delta))) +
  # scale_color_manual(values = c("blue", "red")) +
  # geom_smooth(aes(x = dates, y = delta), color = "black", size = 2) +
  geom_hline(yintercept = 0, color = "black") +
  labs(x = "Date", y = "AOU - AOOU (corrected)", color = "NCP Status") +
  theme_bw()

a <- ggplot(data = aou.df) +
  geom_line(aes(x = dates, y = abs(o2_bio - aou)), color = "black") +
  geom_line(aes(x = dates, y = abs(model.performance)), color = "orange") +
  # scale_color_manual(values = c("blue", "red")) +
  # geom_smooth(aes(x = dates, y = delta), color = "black", size = 2) +
  geom_hline(yintercept = 0, color = "black") +
  labs(x = "Date", y = "Delta") +
  theme_bw()
ggplotly(a)

b <- ggplot(data = aou.df) +
  geom_line(aes(x = dates, y = aou), color = "black") +
  geom_line(aes(x = dates, y = o2_bio), color = "blue") +
  geom_line(aes(x = dates, y = aou.corrected), color = "orange") +
  geom_hline(yintercept = 0, color = "black") +
  labs(x = "Date", y = "NCP") +
  theme_bw()
  
b/a

ggplot(data = aou.df) +
  geom_line(aes(x = dates, y = abs(aou.corrected - aou)), color = "red") +
  # geom_line(aes(x = dates, y = abs(model.performance)), color = "orange") +
  # scale_color_manual(values = c("blue", "red")) +
  # geom_smooth(aes(x = dates, y = delta), color = "black", size = 2) +
  geom_hline(yintercept = 0, color = "black") +
  labs(x = "Date", y = "Delta") +
  theme_bw()


# how much better does model predict AOU?
cor.mean <- mean(abs(aou.df$model.performance[which(is.na(aou.df$model.performance) == FALSE)]))
cor.med <- median(abs(aou.df$model.performance[which(is.na(aou.df$model.performance) == FALSE)]))

summary(abs(aou.df$model.performance[which(is.na(aou.df$model.performance) == FALSE)]))

orig.mean <- mean(abs(aou.df$no.model.performance[which(is.na(aou.df$no.model.performance) == FALSE)]))
orig.med <- median(abs(aou.df$no.model.performance[which(is.na(aou.df$no.model.performance) == FALSE)]))

summary(abs(aou.df$no.model.performance[which(is.na(aou.df$no.model.performance) == FALSE)]))

1-cor.mean/orig.mean
1-cor.med/orig.med




# lmfit <- lm(data = aou.df,
#             aou.corrected ~ sin(2*pi*X/(24)) + cos(2*pi*X/(24)))
# b0 <- coef(lmfit)[1]
# alpha <- coef(lmfit)[2]
# beta <- coef(lmfit)[3]
# r <- sqrt(alpha^2 + beta^2)
# phi <- atan2(beta, alpha)
# 
# summary(lmfit)
# 
# par(mfrow=c(1,2))
# curve(b0 + r * sin(x + phi), 0, 2*pi, lwd=3, col="Gray",
#       main="Overplotted Graphs", xlab="x", ylab="y")
# curve(b0 + alpha * sin(x) + beta * cos(x), lwd=3, lty=3, col="Red", add=TRUE)
# curve(b0 + r * sin(x + phi) - (b0 + alpha * sin(x) + beta * cos(x)), 
#       0, 2*pi, n=257, lwd=3, col="Gray", main="Difference", xlab="x", y="")


# ---- fitting trig function ----

# aou.df$X <- aou.df$X/(24*365.25)

ssp <- spectrum(aou.df$aou.corrected)  
per <- 1/ssp$freq[ssp$spec==max(ssp$spec)]
aou.df$num.date <- as.numeric(difftime(aou.df$dates, min(aou.df$dates), units = "hour"))
reslm <- lm(aou.df$aou.corrected ~ sin(2*pi/per*aou.df$num.date)+cos(2*pi/per*aou.df$num.date))
summary(reslm)

rg <- diff(range(aou.df$aou.corrected))

plot(aou.df$aou.corrected~aou.df$dates, ylim=c(min(aou.df$aou.corrected)-0.1*rg,max(aou.df$aou.corrected)+0.1*rg), pch = 16, cex = 0.5)
lines(fitted(reslm)~aou.df$dates,col=4,lty=2)   # dashed blue line is sin fit

# including 2nd harmonic really improves the fit
reslm2 <- lm(aou.df$aou.corrected ~ sin(2*pi/per*aou.df$dates)+cos(2*pi/per*aou.df$dates)+sin(4*pi/per*aou.df$dates)+cos(4*pi/per*aou.df$dates))
summary(reslm2)
lines(fitted(reslm2)~aou.df$dates,col=3)    # solid green line is periodic with second harmonic





# TRY AGGREGATING IN DIFFERENT WAYS


# rolling avg (BETTER)
try.it <- aou.df[which(is.na(aou.df$aou.cor.rolling.avg) == FALSE),]
ssp <- spectrum(try.it$aou.cor.rolling.avg)  
per <- 1/ssp$freq[ssp$spec==max(ssp$spec)]
try.it$num.date <- as.numeric(difftime(try.it$dates, min(try.it$dates), units = "hour"))
reslm <- lm(try.it$aou.cor.rolling.avg ~ sin(2*pi/per*try.it$num.date))
summary(reslm)

rg <- diff(range(try.it$aou.cor.rolling.avg))

plot(try.it$aou.cor.rolling.avg~try.it$dates, ylim=c(min(try.it$aou.cor.rolling.avg)-0.1*rg,max(try.it$aou.cor.rolling.avg)+0.1*rg), pch = 16, cex = 0.5)
lines(fitted(reslm)~try.it$dates,col="blue")   # dashed blue line is sin fit

per/(24*365.25) # period in years
per/24

# rolling average -- just linear
try.it <- aou.df[which(is.na(aou.df$aou.cor.rolling.avg) == FALSE),]
ssp <- spectrum(try.it$aou.cor.rolling.avg)  
per <- 1/ssp$freq[ssp$spec==max(ssp$spec)]
try.it$num.date <- as.numeric(difftime(try.it$dates, min(try.it$dates), units = "hour"))
reslm <- lm(try.it$aou.cor.rolling.avg ~ try.it$num.date)
summary(reslm)

rg <- diff(range(try.it$aou.cor.rolling.avg))

plot(try.it$aou.cor.rolling.avg~try.it$dates, ylim=c(min(try.it$aou.cor.rolling.avg)-0.1*rg,max(try.it$aou.cor.rolling.avg)+0.1*rg), pch = 16, cex = 0.5)
lines(fitted(reslm)~try.it$dates,col=4,lty=2)   # dashed blue line is sin fit


# rolling avg  -- with linear addition (BEST)
try.it <- aou.df[which(is.na(aou.df$aou.cor.rolling.avg) == FALSE),]
ssp <- spectrum(try.it$aou.cor.rolling.avg)  
per <- 1/ssp$freq[ssp$spec==max(ssp$spec)]
try.it$num.date <- as.numeric(difftime(try.it$dates, min(try.it$dates), units = "hour"))
reslm <- lm(try.it$aou.cor.rolling.avg ~ sin(2*pi/per*try.it$num.date)+cos(2*pi/per*try.it$num.date) + try.it$num.date)
summary(reslm)

rg <- diff(range(try.it$aou.cor.rolling.avg))

plot(try.it$aou.cor.rolling.avg~try.it$dates, ylim=c(min(try.it$aou.cor.rolling.avg)-0.1*rg,max(try.it$aou.cor.rolling.avg)+0.1*rg), pch = 16, cex = 0.5)
lines(fitted(reslm)~try.it$dates,col=4,lty=2)   # dashed blue line is sin fit

per/(24*365.25) # period in years
per/24

0.123+0.2013 # combined trig and linear model fit better than each individually!!









