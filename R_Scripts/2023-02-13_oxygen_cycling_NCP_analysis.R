# Analysis of oxygen cycling and model NCP/AOP predictions
# RJH
# 2023-02-13

# ---- library ----

library(tidyverse)
library(lubridate)
library(zoo)

# ---- read in data ----

setwd("C://Users/haler/Documents/PhD-Bowman/MIMS-miniDOT_O2-Ar_Study/MIMS_4riley/")

model.input.data <- readRDS(file = "2023-02-03_combined.df.select.rds") # this is the dataset used in the model
corrected.full.time.series.data <- readRDS(file = "full_aou_df.rds")
mims.zoo <- readRDS(file = "2023-02-03_mims_zoo.rds")


# ---- functions ----

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
  
  curve(expr = dnorm(x, mean = model.mean, sd = model.sd), col = "blue", from = min((model.mean-(4*model.sd)), (all.mean-(4*all.sd))), to = max((model.mean+(4*model.sd)), (all.mean+(4*all.sd))), ylab = "Frequency", xlab = "[O2]bio / AOU")
  curve(expr = dnorm(x, mean = all.mean, sd = all.sd), col = "red", add = TRUE)
  legend('topleft', legend = c("model", "all"), col = c("blue", "red"), lty = 1)
  
  return(my.summary)
  
}


# ---- format data ----

model.input.data$dates <- parse_date_time(as.character(model.input.data$dates), orders = "Ymd HMS")
corrected.full.time.series.data$dates <- parse_date_time(as.character(corrected.full.time.series.data$dates), orders = "Ymd HMS")

mims.df <- data.frame(mims.zoo)
mims.df$Group.1 <- parse_date_time(as.character(mims.df$Group.1), order = "Ymd HMS")
colnames(mims.df)[1] <- "dates"
mims.df <- mims.df %>% mutate_at(colnames(mims.df)[-1], as.numeric)

# ---- analysis of net trophic status ----

model.dates1 <- min(model.input.data$dates)
model.dates2 <- max(model.input.data$dates)

# MIMS O2_bio
compare.model.to.all.data(model.input.data$o2_bio, mims.df$o2_bio)
t.test(model.input.data$o2_bio, mu = 0, alternative = "less") # statistics show that the o2_bio data is not neutral trophic status

qqnorm(mims.df$o2_bio, pch = 1, frame = FALSE)
qqline(mims.df$o2_bio, col = "steelblue", lwd = 2)


# corrected AOU
# all
t.test(corrected.full.time.series.data$aou.corrected, mu = 0, alternative = "two.sided") # statistics show that the o2_bio data is not neutral trophic status

# only model dates
test.this <- corrected.full.time.series.data[which(corrected.full.time.series.data$dates > model.dates1 & corrected.full.time.series.data$dates < model.dates2),]
t.test(test.this$aou.corrected, mu = 0, alternative = "two.sided") # statistics show that the o2_bio data is not neutral trophic status



# uncorrected AOU
t.test(corrected.full.time.series.data$aou, mu = 0, alternative = "two.sided") # statistics show that the o2_bio data is not neutral trophic status

# only model dates
test.this <- corrected.full.time.series.data[which(corrected.full.time.series.data$dates > model.dates1 & corrected.full.time.series.data$dates < model.dates2),]
t.test(test.this$aou, mu = 0, alternative = "two.sided") # statistics show that the o2_bio data is not neutral trophic status
mean(test.this$aou)



curve(expr = dnorm(x, mean = mean(corrected.full.time.series.data$aou), sd = sd(corrected.full.time.series.data$aou)), col = "blue", from = -100, to = 100)
# abline(v = mean(corrected.full.time.series.data$aou), col = "blue")
curve(expr = dnorm(x, mean = mean(test.this$aou), sd = sd(test.this$aou)), col = "blue", lwd = 4, from = -100, to = 100, add = TRUE)
# abline(v = mean(corrected.full.time.series.data$aou.corrected), col = "red")

curve(expr = dnorm(x, mean = mean(corrected.full.time.series.data$aou.corrected), sd = sd(corrected.full.time.series.data$aou.corrected)), col = "red", from = -100, to = 100, add = TRUE)
# abline(v = mean(corrected.full.time.series.data$aou.corrected), col = "red")
curve(expr = dnorm(x, mean = mean(test.this$aou.corrected), sd = sd(test.this$aou.corrected)), col = "red", lwd = 4, from = -100, to = 100, add = TRUE)
# abline(v = mean(corrected.full.time.series.data$aou.corrected), col = "red")

curve(expr = dnorm(x, mean = mean(model.input.data$o2_bio), sd = sd(model.input.data$o2_bio)), col = "green", from = -100, to = 100, lwd = 4, add = TRUE)
# abline(v = mean(model.input.data$o2_bio), col = "green")


t.test(x = corrected.full.time.series.data$aou, y = corrected.full.time.series.data$aou.corrected)
t.test(x = test.this$aou, y = test.this$aou.corrected)
t.test(x = test.this$aou, y = model.input.data$o2_bio)
t.test(x = test.this$aou.corrected, y = model.input.data$o2_bio)



plot(corrected.full.time.series.data$aou~corrected.full.time.series.data$dates, pch = 16, col = alpha("blue", 0), xlab = "Date", ylab = "O2 Anomaly")
lines(corrected.full.time.series.data$aou~corrected.full.time.series.data$dates, pch = 16, col = alpha("blue", 0.5))
lines(corrected.full.time.series.data$aou.corrected~corrected.full.time.series.data$dates, pch = 16, col = alpha("red", 0.5))
lines(model.input.data$o2_bio~model.input.data$dates, pch = 16, col = alpha("green", 0.5))
abline(h = mean(corrected.full.time.series.data$aou), col = "blue", lwd = 2)
abline(h = mean(corrected.full.time.series.data$aou.corrected), col = "red", lwd = 2)
segments(y0 = mean(model.input.data$o2_bio), y1 = mean(model.input.data$o2_bio), x0 = model.dates1, x1 = model.dates2, col = "green", lwd = 2)
legend("topright", legend = c("AOP", "AOP_cor", "[O2]bio"), col = c("blue", "red", "green"), lwd = 2)

plot(corrected.full.time.series.data$aou~corrected.full.time.series.data$dates, pch = 16, col = alpha("blue", 0), xlab = "Date", ylab = "O2 Anomaly")
lines((corrected.full.time.series.data$aou-corrected.full.time.series.data$aou.corrected)~corrected.full.time.series.data$dates, pch = 16, col = alpha("blue", 0.5))

# ---- model improvement of AOU ----

try.it <- merge(corrected.full.time.series.data, model.input.data, by = c("dates", "aou"))

100 * mean(((try.it$o2_bio - try.it$aou) - (try.it$o2_bio - try.it$aou.corrected)) / (try.it$o2_bio - try.it$aou))











