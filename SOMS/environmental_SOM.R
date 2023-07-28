## This code is for SOM segmentation of environmental DNA collected from the MOSAiC International Arctic Drift Expedition.
#Code written by Emelia Chamberlain, adapted from code written by Jeff Bowman for Bowman et al., 2016, ISMEJ.

#updated Feb 2023

########## Set working directory and load libraries ##########
# Set working directory for output
setwd("C://Users/haler/Documents/PhD-Bowman/MIMS-miniDOT_O2-Ar_Study/SOMS/")

#load libarires
library(kohonen)
library(stringr)
library(gplots)
library(vegan)
library(ggpubr)
library(cowplot)

########## Set up data ##########

load("../16S_sccoos/20230510_sccoos_asv.Rdata")

#running on 16s only (for now) 

## issue with 16s names, fix NEED TO CHECK OUT IN MERGING FILES LATER
x <- str_sub(rownames(fullmeta), end=-8)
y <- paste0(x, "6S.exp.")
fullmeta$libraryname_pro <- y
rm(x,y)

#make rownames 16s library names
rownames(fullmeta) <- fullmeta$libraryname_pro

#reorder to match tally files
fullmeta <- fullmeta[rownames(unique.16s.select),]

## isolate just environmental data
meta.enviro <- subset(fullmeta, type == "FT_DNA" | type == "CTD" | type == "ICE" | type == "DMS ice core" | type == "CTD " | type == "Melt_Pond" | type == "Under_Ice" | type == "Melt_Pond_ICE" | type == "Lead ")
merged.pro <- unique.16s.select[rownames(meta.enviro),]

## Delete ASVs and edges which don't appear in this dataset 
merged.pro <- merged.pro[,which(colSums(merged.pro) != 0)]

## Clean up map files to match 
merged.pro.map <- fullmap.16s[colnames(merged.pro),]
merge.pro.taxa <- taxa.16s[colnames(merged.pro),]

## Transform ASV count data
#Take relative abundance
merged.pro.RA <- merged.pro/rowSums(merged.pro)

# Hellinger transformation, this will take the square root of the relative abundance. 
merged.pro.hell <- decostand(merged.pro, method = "hellinger")

#rename metadata
merged.pro.meta <- meta.enviro


## remove excess data (for now) to clean up environment while working
rm(unique.16s.select, tally.16s, taxa.16s, fullmap.16s)
rm(depth_repository, meta_archive, fullmeta)
rm(fullmap.18s, tally.18s, taxa.18s, unique.18s.select)
rm(euk.meta.colnames, pro.meta.colnames)

## save chloroplast data 
write.csv(chloromap, "~/mosaic/new_2023/Data_QC/DNA/generated_data/chloromap.csv")
write.csv(chlorotally, "~/mosaic/new_2023/Data_QC/DNA/generated_data/chlorotally.csv")

rm(chloromap, chlorotally)

## create matrices for input
#option 1 - hellinger transform
pro.matrix <- as.matrix(merged.pro.hell)
#option 2 - scale function in r. The scale() function subtracts the values of each column by the matching “center” value from the argument. This is also known as data standardization, and it basically involves converting each original value into a z-score.
#If the value is numeric, the scale() method divides the values of each column by the corresponding scale value from the input.
pro.scale <- as.matrix(scale(merged.pro.RA))

########### calculate ESOM and clusters #############

### build ESOM - prokaryotes ONLY (this file)

som.grid <- somgrid(xdim = 10, ydim = 10, topo="hexagonal", toroidal = T)

set.seed(2021)

## hellinger transformed data
som.model.pro <- som(pro.matrix,
                     grid = som.grid,
                     rlen = 100,
                     alpha = c(0.05,0.01),
                     keep.data = TRUE)
som.events.pro <- som.model.pro$codes[[1]]
som.dist.pro <- as.matrix(dist(som.events.pro))

## scaled RA data
som.model.pro.scale <- som(pro.scale,
                           grid = som.grid,
                           rlen = 100,
                           alpha = c(0.05,0.01),
                           keep.data = TRUE)
som.events.pro.scale <- som.model.pro.scale$codes[[1]]
som.dist.pro.scale <- as.matrix(dist(som.events.pro.scale))

## hell transform
plot(som.model.pro, type = 'mapping', pch = 19, palette.name = topo.colors, main = '')
my.data <- som.events.pro
wss.pro <- (nrow(my.data)-1)*sum(apply(my.data,2,var)) 
for (i in 2:24) {
  wss.pro[i] <- sum(kmeans(my.data, centers=i)$withinss)
}

plot(wss.pro)

## scale
plot(som.model.pro.scale, type = 'mapping', pch = 19, palette.name = topo.colors, main = '')
my.data <- som.events.pro.scale
wss.pro.scale <- (nrow(my.data)-1)*sum(apply(my.data,2,var)) 
for (i in 2:24) {
  wss.pro.scale[i] <- sum(kmeans(my.data, centers=i)$withinss)
}

plot(wss.pro.scale)



k1 = 9


som.cluster.pro <- kmeans(som.events.pro, centers = k1)
som.cluster.pro.scale <- kmeans(som.events.pro.scale, centers = k1)


## plots
plot(som.model.pro,
    main = '',
     type = "property",
     property = som.cluster.pro$cluster,
     palette.name = rainbow)
add.cluster.boundaries(som.model.pro, som.cluster.pro$cluster)


plot(som.model.pro.scale,
     main = '',
     type = "property",
     property = som.cluster.pro.scale$cluster,
     palette.name = topo.colors)
add.cluster.boundaries(som.model.pro.scale, som.cluster.pro.scale$cluster)

aweSOMplot(som = som.model.pro.scale, type = "Barplot", data = merged.pro.meta, 
           variables = c("depth_m", "temp", "salinity", "chla", "O2bio"), 
           superclass = som.cluster.pro.scale$cluster)

x <- aweSOMplot(som = som.model.pro.scale, type = "Color", data = merged.pro.meta, 
                variables = c("salinity"), 
                superclass = som.cluster.pro.scale$cluster)

## OK after some messing around, we're going to add a few different options to our metadata file to export for final plotting options, once we look at the transitions hopefully things will make more sense?


##scaled
som.cluster9.pro.scale <- kmeans(som.events.pro.scale, centers = 9)
merged.pro.meta$pro_mode_scale_k9 <- som.cluster9.pro.scale$cluster[som.model.pro.scale$unit.classif]

som.cluster10.pro.scale <- kmeans(som.events.pro.scale, centers = 10)
merged.pro.meta$pro_mode_scale_k10 <- som.cluster10.pro.scale$cluster[som.model.pro.scale$unit.classif]

som.cluster11.pro.scale <- kmeans(som.events.pro.scale, centers = 11)
merged.pro.meta$pro_mode_scale_k11 <- som.cluster11.pro.scale$cluster[som.model.pro.scale$unit.classif]

som.cluster12.pro.scale <- kmeans(som.events.pro.scale, centers = 12)
merged.pro.meta$pro_mode_scale_k12 <- som.cluster11.pro.scale$cluster[som.model.pro.scale$unit.classif]

##hell
som.cluster9.pro <- kmeans(som.events.pro, centers = 9)
merged.pro.meta$pro_mode_hell_k9 <- som.cluster9.pro$cluster[som.model.pro$unit.classif]

som.cluster10.pro <- kmeans(som.events.pro, centers = 10)
merged.pro.meta$pro_mode_hell_k10 <- som.cluster10.pro$cluster[som.model.pro$unit.classif]

som.cluster11.pro <- kmeans(som.events.pro, centers = 11)
merged.pro.meta$pro_mode_hell_k11 <- som.cluster11.pro$cluster[som.model.pro$unit.classif]

########## testing plots ##########

##fix date issue
nums <- c(177,178,179,180,181,182,183)
dates <- as.Date(c("2020-03-28", "2020-03-28", "2020-03-28", "2020-03-28", "2020-03-28", "2020-03-28", "2020-03-28"))

##slim type for plotting
merged.pro.meta[which(merged.pro.meta$type == "ICE"), 8] <- "Ice"
merged.pro.meta[which(merged.pro.meta$type == "Melt_Pond_ICE"), 8] <- "Ice"
merged.pro.meta[which(merged.pro.meta$type == "DMS ice core"), 8] <- "Ice"
merged.pro.meta[which(merged.pro.meta$type == "CTD "), 8] <- "CTD"
unique(merged.pro.meta$type)

##add fake depth for ice 

merged.pro.meta[which(merged.pro.meta$type == "Ice"), 54] <- -5


for (i in 1:length(nums)) {
  merged.pro.meta[nums[i],47] <- dates[i]
}

aweSOMplot(som = som.model.pro, type = "Color", data = merged.pro.meta, variables = "salinity")

z <- rownames(merged.pro.meta)
#as.character(pro_mode_hell_k11)
ggplot(merged.pro.meta, aes(x=date, y = depth_m, color= salinity)) + geom_point() + ylim(100,-8)

## move forward with pro mode hell 11? 

plot(som.model.pro.scale,
     main = '',
     type = "property",
     property = som.cluster11.pro$cluster,
     palette.name = topo.colors)
add.cluster.boundaries(som.model.pro, som.cluster11.pro$cluster)

## save data and move to final plots

ggplot(merged.pro.meta, aes(x=as.character(pro_mode_scale_k10), y = depth_m)) + geom_boxplot() + theme_light() + xlab("Type") + ylab("Temperature") + theme(legend.position = "top",  legend.title = element_blank(), legend.background = element_rect(size = 0.5, linetype = "solid", colour = "grey"), axis.title = element_text(size=16), axis.text = element_text(size = 14)) 
#add dots and distinguish between sample types
temp <- temp + geom_beeswarm(aes(color = Type), size = 0.5) 
#add stats
temp + stat_compare_means(aes(group = melt_type), label = "p.signif", show.legend = FALSE, size = 5)
