# MIMS predictor ordination
# 2022-10-03
# RJH

# ---- library ----

library(tidyverse)
library(vegan)
library(lubridate)
library(plotly)


# ---- import data ---

setwd("C://Users/haler/Documents/PhD-Bowman/MIMS-miniDOT_O2-Ar_Study/MIMS_4riley/")

full.aou.df <- readRDS("full_aou_df.rds")
full.predictors.df <- readRDS("full_predictors_df.rds")
mims.hourly <- readRDS("mims_hourly.rds")

# ---- create dissimilarity matrix ----

my.dates.aou <- full.predictors.df[,c("aou", "dates")]
full.predictors.df <- full.predictors.df[,c("temperature.sccoos.zoo", "salinity", "pressure", "WSPD", "GST", "WVHT", "DPD")]
rownames(full.predictors.df) <- as.character(my.dates.aou$dates)
predictors.matrix <- as.matrix(full.predictors.df)

predictors.matrix <- scale(predictors.matrix)


# ---- NMDS ----
distmat <- vegdist(predictors.matrix, method = "bray", binary = F, diag = F, upper = F, na.rm = T)

NMS <- metaMDS(distmat, distance = "bray")
goodness(NMS)
stressplot(NMS)

plot(NMS, type = "t")

data.scores <- as.data.frame((NMS$points))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$Sample.Name <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores


# match PCA so I don't have to redo it
df1 <- data.scores
df1 <- merge(df1, save.for.later, by = "Sample.Name", all = FALSE)

colnames(df1)[c(2,3)] <- c("PC1", "PC2")


# ---- PCA ----
my.pca <- rda(predictors.matrix)
smry <- summary(my.pca)
df1  <- data.frame(smry$sites[,1:2])       # PC1 and PC2
df2  <- data.frame(smry$species[,1:2])

df1$Date <- rownames(df1)
df1$Date <- parse_date_time(df1$Date, orders = "Ymd HMS")

colnames(full.aou.df)[3] <- "Date"
full.aou.df$Date <- parse_date_time(full.aou.df$Date, orders = "Ymd HMS")
df1 <- merge(df1, full.aou.df, by = "Date")
full.predictors.df$Date <- parse_date_time(rownames(full.predictors.df), orders = "Ymd HMS")
df1 <- merge(df1, full.predictors.df, by = "Date")

colnames(mims.hourly)[1] <- "Date"
mims.hourly$Date <- parse_date_time(mims.hourly$Date, orders = "Ymd HMS")
o2_bio_hourly <- mims.hourly[,c("Date", "o2_bio")]
df1 <- merge(df1, o2_bio_hourly, by = "Date", all = TRUE)

df1$delta <- df1$aou - df1$aou.corrected

df1$decile <- ntile(abs(df1$delta), 10)


# ---- plot----

# AOU-AOUcor
ggplot(data = df1, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = (abs(delta))), size = 3, alpha = 0.2) +
  # geom_point(aes(color = log(abs(delta))), size = 3, alpha = 0.2) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  # scale_color_continuous(type = "viridis", name = "AOU - AOU(corrected)") + 
  scale_color_viridis_c(option = "inferno", name = "log(abs(AOU - AOU(corrected)))") +
  # coord_fixed() +
  geom_segment(data=df2, aes(x=0, xend=PC1/10, y=0, yend=PC2/10),
               color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2,
            aes(x=PC1/10,y=PC2/10,label=rownames(df2),
                ),
            color="red", size=4) +
  theme_bw()

# AOU-AOUcor
ggplot(data = df1, aes(x = PC1, y = PC2)) +
  # geom_point(aes(color = (abs(delta))), size = 3, alpha = 0.2) +
  geom_point(aes(color = log(delta)), size = 3, alpha = 0.2) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  # scale_color_continuous(type = "viridis", name = "AOU - AOU(corrected)") + 
  scale_color_viridis_c(option = "inferno", name = "log(AOU - AOU(corrected))") +
  # coord_fixed() +
  geom_segment(data=df2, aes(x=0, xend=PC1/10, y=0, yend=PC2/10),
               color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2,
            aes(x=PC1/10,y=PC2/10,label=rownames(df2),
            ),
            color="red", size=4) +
  theme_bw()

# AOU-AOUcor
ggplot(data = df1, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = (abs(delta))), size = 3, alpha = 0.2) +
  # geom_point(aes(color = log(abs(delta))), size = 3, alpha = 0.2) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  # scale_color_continuous(type = "viridis", name = "AOU - AOU(corrected)") + 
  scale_color_viridis_c(option = "inferno", name = "abs(AOU - AOU(corrected))") +
  # coord_fixed() +
  geom_segment(data=df2, aes(x=0, xend=PC1/10, y=0, yend=PC2/10),
               color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2,
            aes(x=PC1/10,y=PC2/10,label=rownames(df2),
            ),
            color="red", size=4) +
  theme_bw()

# AOU-AOUcor
ggplot(data = df1, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = (delta)), size = 3, alpha = 0.2) +
  # geom_point(aes(color = log(abs(delta))), size = 3, alpha = 0.2) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  # scale_color_continuous(type = "viridis", name = "AOU - AOU(corrected)") + 
  scale_color_viridis_c(option = "inferno", name = "AOU - AOU(corrected)") +
  # coord_fixed() +
  geom_segment(data=df2, aes(x=0, xend=PC1/10, y=0, yend=PC2/10),
               color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2,
            aes(x=PC1/10,y=PC2/10,label=rownames(df2),
            ),
            color="red", size=4) +
  theme_bw()
# extreme delta (difference between AOU and AOU_corrected) driven by temperature
# sample points with high temperatures are corrected more in the model
# measured AOU values are more incorrect in sample points with high temperatures 

ggplot(data = df1) +
  geom_point(aes(x = temperature.sccoos.zoo, y = delta, color = (abs(delta))), size = 3, alpha = 0.2) +
  scale_color_viridis_c(option = "inferno", name = "abs(AOU - AOU(corrected))") +
  # coord_fixed() +
  theme_bw()
# nope nevermind, don\t see that here


df1.subset <- df1[df1$decile == 10,]

ggplot(data = df1.subset, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = abs(delta)), size = 3, alpha = 0.2) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  # scale_color_continuous(type = "viridis", name = "AOU - AOU(corrected)") + 
  scale_color_viridis_c(option = "inferno", name = "abs(AOU - AOU(corrected))") +
  # coord_fixed() +
  geom_segment(data=df2, aes(x=0, xend=PC1/10, y=0, yend=PC2/10),
               color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2,
            aes(x=PC1/10,y=PC2/10,label=rownames(df2),
            ),
            color="red", size=4) +
  theme_bw()
# less evident just looking at the extreme delta values

ggplot(data = df1.subset, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = delta), size = 3, alpha = 0.4) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  # scale_color_continuous(type = "viridis", name = "AOU - AOU(corrected)") + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  # coord_fixed() +
  geom_segment(data=df2, aes(x=0, xend=PC1/10, y=0, yend=PC2/10),
               color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2,
            aes(x=PC1/10,y=PC2/10,label=rownames(df2),
            ),
            color="red", size=4) +
  theme_bw()
# higher environmental drivers (temp, DPD, GST/WSPD) important in samples where measured AOU is lower than AOU corrected by model
# is this the right interpretation???



# colored by environmental variables 
# duh of course these are correlated... at least they look pretty
ggplot(data = df1, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = temperature.sccoos.zoo), size = 3, alpha = 0.2) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  # scale_color_continuous(type = "viridis", name = "AOU - AOU(corrected)") + 
  scale_color_viridis_c(option = "inferno", name = "Temperature") +
  # coord_fixed() +
  geom_segment(data=df2, aes(x=0, xend=PC1/10, y=0, yend=PC2/10),
               color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2,
            aes(x=PC1/10,y=PC2/10,label=rownames(df2),
            ),
            color="red", size=4) +
  theme_bw()

ggplot(data = df1, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = DPD), size = 3, alpha = 0.2) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  # scale_color_continuous(type = "viridis", name = "AOU - AOU(corrected)") + 
  scale_color_viridis_c(option = "inferno", name = "DPD") +
  # coord_fixed() +
  geom_segment(data=df2, aes(x=0, xend=PC1/10, y=0, yend=PC2/10),
               color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2,
            aes(x=PC1/10,y=PC2/10,label=rownames(df2),
            ),
            color="red", size=4) +
  theme_bw()

ggplot(data = df1, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = GST), size = 3, alpha = 0.2) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  # scale_color_continuous(type = "viridis", name = "AOU - AOU(corrected)") + 
  scale_color_viridis_c(option = "inferno", name = "GST") +
  # coord_fixed() +
  geom_segment(data=df2, aes(x=0, xend=PC1/10, y=0, yend=PC2/10),
               color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2,
            aes(x=PC1/10,y=PC2/10,label=rownames(df2),
            ),
            color="red", size=4) +
  theme_bw()




# just negative model AOU corrections: measured AOU higher than model (model decreased AOU to match [O2]bio)
df1.subset <- df1[which(df1$delta > 0),]

ggplot(data = df1.subset, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = delta), size = 3, alpha = 0.2) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  # scale_color_continuous(type = "viridis", name = "AOU - AOU(corrected)") + 
  scale_color_viridis_c(option = "inferno", name = "abs(AOU - AOU(corrected))") +
  # coord_fixed() +
  geom_segment(data=df2, aes(x=0, xend=PC1/10, y=0, yend=PC2/10),
               color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2,
            aes(x=PC1/10,y=PC2/10,label=rownames(df2),
            ),
            color="red", size=4) +
  theme_bw()
# high temp samples were more changed by model (decreasing AOU)

ggplot(data = df1.subset) +
  geom_point(aes(x = temperature.sccoos.zoo, y = delta, color = (abs(delta))), size = 3, alpha = 0.2) +
  scale_color_viridis_c(option = "inferno", name = "AOU - AOU(corrected)") +
  # coord_fixed() +
  theme_bw()
# oooooo here's something maybe, or is this just an artifact of looking at part of the data
# 


# just positive model AOU corrections 
df1.subset <- df1[which(df1$delta < 0),]

ggplot(data = df1.subset, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = -delta), size = 3, alpha = 0.2) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  # scale_color_continuous(type = "viridis", name = "AOU - AOU(corrected)") + 
  scale_color_viridis_c(option = "inferno", name = "abs(AOU - AOU(corrected))") +
  # coord_fixed() +
  geom_segment(data=df2, aes(x=0, xend=PC1/10, y=0, yend=PC2/10),
               color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2,
            aes(x=PC1/10,y=PC2/10,label=rownames(df2),
            ),
            color="red", size=4) +
  theme_bw()
# hmm don't really see the same thing here

ggplot(data = df1.subset) +
  geom_point(aes(x = temperature.sccoos.zoo, y = delta, color = (abs(delta))), size = 3, alpha = 0.2) +
  scale_color_viridis_c(option = "inferno", name = "AOU - AOU(corrected)") +
  # coord_fixed() +
  theme_bw()
# maybe something here? why are positive and negative values behaving differently?




df1.subset <- df1[df1$decile == 10,]

ggplot(data = df1.subset, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = abs(aou)), size = 3, alpha = 0.2) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  # scale_color_continuous(type = "viridis", name = "AOU - AOU(corrected)") + 
  scale_color_viridis_c(option = "inferno", name = "AOU") +
  # coord_fixed() +
  geom_segment(data=df2, aes(x=0, xend=PC1/10, y=0, yend=PC2/10),
               color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2,
            aes(x=PC1/10,y=PC2/10,label=rownames(df2),
            ),
            color="red", size=4) +
  theme_bw()
# not a very clear correlation. High temp samples do not correlate with extreme AOU (measured)

ggplot(data = df1, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = abs(aou)), size = 3, alpha = 0.2) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  # scale_color_continuous(type = "viridis", name = "AOU - AOU(corrected)") + 
  scale_color_viridis_c(option = "inferno", name = "AOU") +
  # coord_fixed() +
  geom_segment(data=df2, aes(x=0, xend=PC1/10, y=0, yend=PC2/10),
               color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2,
            aes(x=PC1/10,y=PC2/10,label=rownames(df2),
            ),
            color="red", size=4) +
  theme_bw()
# yep still no correlation popping out to me

ggplot(data = df1, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = abs(aou.corrected)), size = 3, alpha = 0.2) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  # scale_color_continuous(type = "viridis", name = "AOU - AOU(corrected)") + 
  scale_color_viridis_c(option = "inferno", name = "AOU(corrected)") +
  # coord_fixed() +
  geom_segment(data=df2, aes(x=0, xend=PC1/10, y=0, yend=PC2/10),
               color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2,
            aes(x=PC1/10,y=PC2/10,label=rownames(df2),
            ),
            color="red", size=4) +
  theme_bw()
# extreme AOU_corrected driven by temperature

ggplot(data = df1.subset) +
  geom_point(aes(x = temperature.sccoos.zoo, y = aou.corrected, color = (abs(aou.corrected))), size = 3, alpha = 0.2) +
  scale_color_viridis_c(option = "inferno", name = "AOU - AOU(corrected)") +
  # coord_fixed() +
  theme_bw()
# okay I can maybe see the triangle shape


ggplot(data = df1, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = abs(aou.corrected)), size = 3, alpha = 0.2) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  # scale_color_continuous(type = "viridis", name = "AOU - AOU(corrected)") + 
  scale_color_viridis_c(option = "inferno", name = "AOU(corrected)") +
  # coord_fixed() +
  geom_segment(data=df2, aes(x=0, xend=PC1/10, y=0, yend=PC2/10),
               color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2,
            aes(x=PC1/10,y=PC2/10,label=rownames(df2),
            ),
            color="red", size=4) +
  theme_bw()
# and yep, still see a correlation here with high temps and corrected AOU



# difference between [o2]bio and aou and aou_corrected

df1$o2bio.AOU <- df1$o2_bio - df1$aou
df1$o2bio.AOU.cor <- df1$o2_bio - df1$aou.corrected

df1 <- df1[which(is.na(df1$o2bio.AOU) == FALSE),]
df1 <- df1[which(is.na(df1$o2bio.AOU.cor) == FALSE),]

ggplot(data = df1, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = log(abs(o2bio.AOU))), size = 3, alpha = 0.2) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  # scale_color_continuous(type = "viridis", name = "AOU - AOU(corrected)") + 
  scale_color_viridis_c(option = "inferno", name = "[O2]bio - AOU") +
  # coord_fixed() +
  geom_segment(data=df2, aes(x=0, xend=PC1/10, y=0, yend=PC2/10),
               color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2,
            aes(x=PC1/10,y=PC2/10,label=rownames(df2),
            ),
            color="red", size=4) +
  theme_bw()
# maybe a correlation between samples with high temp and high difference between [o2]bio and aou
# is this what our model predictor importance told us already?

ggplot(data = df1, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = log(abs(o2bio.AOU.cor))), size = 3, alpha = 0.2) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  # scale_color_continuous(type = "viridis", name = "AOU - AOU(corrected)") + 
  scale_color_viridis_c(option = "inferno", name = "[O2]bio - AOU(corrected)") +
  # coord_fixed() +
  geom_segment(data=df2, aes(x=0, xend=PC1/10, y=0, yend=PC2/10),
               color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2,
            aes(x=PC1/10,y=PC2/10,label=rownames(df2),
            ),
            color="red", size=4) +
  theme_bw()
# no correlation here- is this a good thing?







# ---- o2 bio vs aou vs aou_cor ----

df1$perc.changed <- 100 * (((df1$o2_bio - df1$aou) - (df1$o2_bio - df1$aou.corrected))/(df1$o2_bio - df1$aou))
df1$amt.changed <- ((df1$o2_bio - df1$aou) - (df1$o2_bio - df1$aou.corrected))

new.df <- df1[which(is.na(df1$amt.changed) == FALSE),]

ggplot(data = new.df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = amt.changed), size = 3, alpha = 0.2) +
  # geom_point(aes(color = log(abs(delta))), size = 3, alpha = 0.2) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  # scale_color_continuous(type = "viridis", name = "AOU - AOU(corrected)") + 
  scale_color_viridis_c(option = "inferno", name = "AOU - AOU(corrected)") +
  # coord_fixed() +
  geom_segment(data=df2, aes(x=0, xend=PC1/10, y=0, yend=PC2/10),
               color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2,
            aes(x=PC1/10,y=PC2/10,label=rownames(df2),
            ),
            color="red", size=4) +
  theme_bw()

summary(lm(df1$amt.changed~df1$temperature.sccoos.zoo))
plot(df1$amt.changed~df1$temperature.sccoos.zoo)

summary(lm(df1$amt.changed~df1$salinity))
plot(df1$amt.changed~df1$salinity)

summary(lm(df1$amt.changed~df1$pressure))
plot(df1$amt.changed~df1$pressure)

summary(lm(df1$amt.changed~df1$WSPD))
plot(df1$amt.changed~df1$WSPD)


summary(lm(df1$temperature.sccoos.zoo~df1$PC1))
summary(lm(df1$temperature.sccoos.zoo~df1$PC2))

summary(lm(df1$salinity~df1$PC1))
summary(lm(df1$salinity~df1$PC2))

summary(lm(df1$WSPD~df1$PC1))
summary(lm(df1$WSPD~df1$PC2))





