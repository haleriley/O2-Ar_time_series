




setwd("~/mosaic/new_2023/all_SOM")
load("SOM_output_202302.RData")


##### final plots for pro k #####
col = brewer.pal(11, "Paired")
pdf('figs/16S_som_property_map_clusters_k11_hell.pdf', width = 5, height = 5)
plot(som.model.pro,
     main = '',
     type = "property",
     property = som.cluster11.pro$cluster,
     palette.name = rainbow)
add.cluster.boundaries(som.model.pro, som.cluster11.pro$cluster)
dev.off()

pdf('figs/16Ssom_mapping_map_k11_hell_clusters.pdf', width = 5, height = 5)
plot(som.model.pro, type = 'counts', pch = 19, palette.name = rainbow, main = '')
add.cluster.boundaries(som.model.pro, som.cluster11.pro$cluster)
dev.off()

pdf('figs/16Ssom_mapping_map_k11_hell.pdf', width = 5, height = 5)
plot(som.model.pro, type = 'mapping', pch = 19, palette.name = rainbow, main = '')
dev.off()

pdf('figs/16S_wss_hell_k11.pdf',
    width = 5,
    height = 5)
plot(wss.pro,
     pch = 19,
     ylab = 'Within-clusters sum of squares',
     xlab = 'K')
abline(h=wss.pro.scale[11])
abline(v=11)
dev.off()

## screenshotted interactive plots
aweSOMplot(som = som.model.pro, type = "Color", data = merged.pro.meta, 
           variables = c("salinity"), 
           superclass = som.cluster11.pro$cluster)

aweSOMplot(som = som.model.pro, type = "Color", data = merged.pro.meta, 
           variables = c("temp"), 
           superclass = som.cluster11.pro$cluster)

aweSOMplot(som = som.model.pro, type = "Color", data = merged.pro.meta, 
           variables = c("O2bio"), 
           superclass = som.cluster11.pro$cluster)

aweSOMplot(som = som.model.pro, type = "Color", data = merged.pro.meta, 
           variables = c("bac_genome_size.mean"), 
           superclass = som.cluster11.pro$cluster)

aweSOMplot(som = som.model.pro, type = "Color", data = merged.pro.meta, 
           variables = c("bac_gRodon.d.mean"), 
           superclass = som.cluster11.pro$cluster)

aweSOMplot(som = som.model.pro, type = "Color", data = merged.pro.meta, 
           variables = c("bac_n16S.mean"), 
           superclass = som.cluster11.pro$cluster)

aweSOMplot(som = som.model.pro, type = "Color", data = merged.pro.meta, 
           variables = c("bac_GC.mean"), 
           superclass = som.cluster11.pro$cluster)

aweSOMplot(som = som.model.pro, type = "Color", data = merged.pro.meta, 
           variables = c("alpha"), 
           superclass = som.cluster11.pro$cluster)


######### section plots ############

cols1 <- c("white", "darkorange", "gold", "lightgreen", "darkgreen", "cyan", "lightblue", "blue", "purple", "violet", "maroon")
cols2 <- c("red", "darkorange", "gold", "lightgreen", "darkgreen", "cyan", "lightblue", "blue", "purple", "violet", "maroon")


ggplot(merged.pro.meta, aes(x=date, y = depth_m, color= as.character(pro_mode_hell_k11))) + geom_point()  + scale_color_manual(values = cols) + theme_classic() # full map


pdf('figs/16Ssom_k11_hell_55m.pdf', width = 10, height = 6)
ggplot(merged.pro.meta, aes(x=date, y = depth_m, color= as.character(pro_mode_hell_k11))) + geom_point(size = 4)  + scale_color_manual(values = cols1) + theme_classic() + ylim(55,-8) + theme(legend.position = "bottom", legend.text = element_text(color = 16), legend.title = element_text(size = 16), axis.text = element_text(size = 16), axis.title = element_text(size = 16)) + labs(color = "Taxonomic Mode", y = "Depth (m)", x = "Date")
dev.off()


pdf('figs/16Ssom_k11_hell_latitude.pdf', width = 10, height = 6)
ggplot(merged.pro.meta, aes(x=date, y = lat, color= as.character(pro_mode_hell_k11))) + geom_point(size = 4)  + scale_color_manual(values = cols2) + theme_classic() + theme(legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10), axis.title = element_text(size = 16)) + labs(color = "Taxonomic Mode", y = "Latitude", x = "Date")
dev.off()


########## nmds #########

## Merged NMDS plot of relative abundances
set.seed(151)
bac.mds <- metaMDS(merged.pro.hell, k = 3) #ordinate, stress = 0.08621695
bac.mds.samples <- bac.mds$points

#prelim plot
pdf("trialplot.pdf")
plot(bac.mds.samples[,1], bac.mds.samples[,2],
     ylab = 'Dim 2',
     xlab = 'Dim 1')
dev.off()
    

pro.mds.samples.df <- as.data.frame(bac.mds.samples)
pro.mds.samples.df$date <- merged.pro.meta$date
pro.mds.samples.df$lat <- merged.pro.meta$lat
pro.mds.samples.df$Type <- merged.pro.meta$type
pro.mds.samples.df$Depth <- merged.pro.meta$depth_m
pro.mds.samples.df$temp <- merged.pro.meta$temp
pro.mds.samples.df$salin <- merged.pro.meta$salinity
pro.mds.samples.df$mode <- merged.pro.meta$pro_mode_hell_k11 


pdf("figs/nmdsmode.pdf")
ggplot(pro.mds.samples.df, aes(x= MDS1, y= MDS2)) + geom_point(aes(MDS1, MDS2, color = as.character(mode), size = date)) + coord_fixed() + theme_classic() + labs(color = "Taxonomic Mode", size = "Date") + theme(legend.text = element_text(size = 16), legend.title = element_text(size = 16), axis.text = element_text(size = 16), axis.title = element_text(size = 16))
dev.off()

pdf("figs/nmdsmodedepth.pdf")
ggplot(pro.mds.samples.df, aes(x= MDS1, y= MDS2)) + geom_point(aes(MDS1, MDS2, color = as.character(mode), size = Depth)) + coord_fixed() + theme_classic() + labs(color = "Taxonomic Mode", size = "Depth (m)") + theme(legend.text = element_text(size = 16), legend.title = element_text(size = 16), axis.text = element_text(size = 16), axis.title = element_text(size = 16)) + scale_size("Depth (m)",breaks=c(-5,0,100,500,4000),labels=c("ice",0,100,500,4000))
dev.off()


pdf("figs/nmdstype.pdf")
ggplot(pro.mds.samples.df, aes(x= MDS1, y= MDS2)) + geom_point(aes(MDS1, MDS2, color = Type, size = date)) + coord_fixed() + theme_classic() + labs(color = "Sample Type", size = "Date") + theme(legend.text = element_text(size = 16), legend.title = element_text(size = 16), axis.text = element_text(size = 16), axis.title = element_text(size = 16))
dev.off()

pdf("figs/nmdstypedepth.pdf")
ggplot(pro.mds.samples.df, aes(x= MDS1, y= MDS2)) + geom_point(aes(MDS1, MDS2, color = Type, size = Depth)) + coord_fixed() + theme_classic() + labs(color = "Sample Type", size = "Depth (m)") + theme(legend.text = element_text(size = 16), legend.title = element_text(size = 16), axis.text = element_text(size = 16), axis.title = element_text(size = 16)) + scale_size("Depth (m)",breaks=c(-5,0,100,500,4000),labels=c("ice",0,100,500,4000))
dev.off()


## which taxa drive variation in MDS plot?
pro.mds.species <- bac.mds$species
pro.mds.species <- as.data.frame(pro.mds.species)

#dimension 1
pro.dim1 <- pro.mds.species[order(-abs(pro.mds.species$MDS1)),]
pro.target.asv.dim1 <- rownames(pro.dim1)
pro.target.asv.dim1 <- pro.target.asv.dim1[1:10]
pro.target.clade.dim1 <- merged.pro.map[pro.target.asv.dim1, "paprica_taxon"]

#dimension 2
pro.dim2 <- pro.mds.species[order(-abs(pro.mds.species$MDS2)),]
pro.target.asv.dim2 <- rownames(pro.dim2)
pro.target.asv.dim2 <- pro.target.asv.dim2[1:10]
pro.target.clade.dim2 <- merged.pro.map[pro.target.asv.dim2, "paprica_taxon"]


######### boxplots ##########

library(ggpubr)

pdf("figs/boxplots/temp.pdf", width = 10, height = 5)
ggplot(merged.pro.meta, aes(x=as.character(pro_mode_hell_k11), y = temp, color = type)) + geom_boxplot() + theme_light() + xlab("Taxonomic Mode") + ylab("Temperature (C)") + theme(legend.position = c(0.9, 0.8),  legend.title = element_blank(), legend.background = element_rect(size = 0.5, linetype = "solid", colour = "grey"), axis.title = element_text(size=16), axis.text = element_text(size = 14)) + stat_compare_means() 
dev.off()


pdf("figs/boxplots/salinity.pdf", width = 10, height = 5)
ggplot(merged.pro.meta, aes(x=as.character(pro_mode_hell_k11), y = salinity, color = type)) + geom_boxplot() + theme_light() + xlab("Taxonomic Mode") + ylab("Salinity") + theme(legend.position = c(0.9, 0.2),  legend.title = element_blank(), legend.background = element_rect(size = 0.5, linetype = "solid", colour = "grey"), axis.title = element_text(size=16), axis.text = element_text(size = 14)) + stat_compare_means() 
dev.off()

pdf("figs/boxplots/o2bio.pdf", width = 10, height = 5)
ggplot(merged.pro.meta, aes(x=as.character(pro_mode_hell_k11), y = O2_Ar*100, color = type)) + geom_boxplot() + theme_light() + xlab("Taxonomic Mode") + ylab("Biologial Oxygen Anomaly (%)") + theme(legend.position = c(0.9, 0.2),  legend.title = element_blank(), legend.background = element_rect(size = 0.5, linetype = "solid", colour = "grey"), axis.title = element_text(size=16), axis.text = element_text(size = 14)) + stat_compare_means() + ylim(-15, 20)
dev.off()

pdf("figs/boxplots/chla.pdf", width = 10, height = 5)
ggplot(merged.pro.meta, aes(x=as.character(pro_mode_hell_k11), y = chla, color = type)) + geom_boxplot() + theme_light() + xlab("Taxonomic Mode") + ylab("Chlorophyll-a") + theme(legend.position = c(0.1, 0.9),  legend.title = element_blank(), legend.background = element_rect(size = 0.5, linetype = "solid", colour = "grey"), axis.title = element_text(size=16), axis.text = element_text(size = 14)) + stat_compare_means() 
dev.off()

pdf("figs/boxplots/methox.pdf", width = 10, height = 5)
ggplot(merged.pro.meta, aes(x=as.character(pro_mode_hell_k11), y = methane_ox, color = type)) + geom_boxplot() + theme_light() + xlab("Taxonomic Mode") + ylab("methane oxidation potential") + theme(legend.position = c(0.1, 0.9),  legend.title = element_blank(), legend.background = element_rect(size = 0.5, linetype = "solid", colour = "grey"), axis.title = element_text(size=16), axis.text = element_text(size = 14)) + stat_compare_means() +ylim(-0.002,0.015)
dev.off()

merged.pro.meta$pro_mode_hell_k11 <- as.character(merged.pro.meta$pro_mode_hell_k11)
#bu <- merged.pro.meta$pro_mode_hell_k11
kt <-  merged.pro.meta %>% kruskal_test(bac_genome_size.mean ~ pro_mode_hell_k11)
pdf("figs/boxplots/bacgnsz.pdf", width = 10, height = 5)
ggplot(merged.pro.meta, aes(x= pro_mode_hell_k11, y = bac_genome_size.mean, color = pro_mode_hell_k11)) + geom_boxplot() + theme_light() + xlab("Taxonomic Mode") + ylab("Estimated bacteria genome size") + theme(legend.position = "none",  legend.title = element_blank(), axis.title = element_text(size=16), axis.text = element_text(size = 14)) + scale_color_manual(values = cols2) + labs(subtitle = get_test_label(kt, detailed = TRUE))
dev.off()

kt <-  merged.pro.meta %>% kruskal_test(bac_GC.mean ~ pro_mode_hell_k11)
pdf("figs/boxplots/bacGC.pdf", width = 10, height = 5)
ggplot(merged.pro.meta, aes(x= pro_mode_hell_k11, y = bac_GC.mean, color = pro_mode_hell_k11)) + geom_boxplot() + theme_light() + xlab("Taxonomic Mode") + ylab("Estimated bacteria GC content") + theme(legend.position = "none",  legend.title = element_blank(), axis.title = element_text(size=16), axis.text = element_text(size = 14)) + scale_color_manual(values = cols2) + labs(subtitle = get_test_label(kt, detailed = TRUE))
dev.off()

kt <-  merged.pro.meta %>% kruskal_test(bac_GC.mean ~ pro_mode_hell_k11)
pdf("figs/boxplots/bacGC.pdf", width = 10, height = 5)
ggplot(merged.pro.meta, aes(x= pro_mode_hell_k11, y = bac_GC.mean, color = pro_mode_hell_k11)) + geom_boxplot() + theme_light() + xlab("Taxonomic Mode") + ylab("Estimated bacteria GC content") + theme(legend.position = "none",  legend.title = element_blank(), axis.title = element_text(size=16), axis.text = element_text(size = 14)) + scale_color_manual(values = cols2) + labs(subtitle = get_test_label(kt, detailed = TRUE))
dev.off()

kt <-  merged.pro.meta %>% kruskal_test(bac_gRodon.d.mean ~ pro_mode_hell_k11)
pdf("figs/boxplots/bacgrodon.pdf", width = 10, height = 5)
ggplot(merged.pro.meta, aes(x= pro_mode_hell_k11, y = bac_gRodon.d.mean, color = pro_mode_hell_k11)) + geom_boxplot() + theme_light() + xlab("Taxonomic Mode") + ylab("Estimated minimal doubling time (hr)") + theme(legend.position = "none",  legend.title = element_blank(), axis.title = element_text(size=16), axis.text = element_text(size = 14)) + scale_color_manual(values = cols2) + labs(subtitle = get_test_label(kt, detailed = TRUE))
dev.off()

kt <-  merged.pro.meta %>% kruskal_test(bac_n16S.mean~ pro_mode_hell_k11)
pdf("figs/boxplots/bac16s.pdf", width = 10, height = 5)
ggplot(merged.pro.meta, aes(x= pro_mode_hell_k11, y = bac_n16S.mean, color = pro_mode_hell_k11)) + geom_boxplot() + theme_light() + xlab("Taxonomic Mode") + ylab("Estimated 16S copy number") + theme(legend.position = "none",  legend.title = element_blank(), axis.title = element_text(size=16), axis.text = element_text(size = 14)) + scale_color_manual(values = cols2) + labs(subtitle = get_test_label(kt, detailed = TRUE))
dev.off()



########## plotting ASV tallies ##########

merged.pro.meta %>% count(type)
z <- merged.pro.meta %>% count(pro_mode_hell_k11, type)

pdf("figs/types.pdf", width = 5, height = 4)
ggplot(data=z, aes(x=n, y=pro_mode_hell_k11, fill=type)) + geom_bar(stat="identity") + theme_classic() + theme(legend.position = c(0.9,0.8), legend.title = element_blank()) + xlab("Taxonomic Mode") + ylab("Number of samples")
dev.off()

## percentage archea vs. bacteria
merged.pro.map %>% count(paprica_domain)

meta$alpha <- diversity(tally_corr[,rownames(meta)],
                        MARGIN = 2,
                        index = "invsimpson")

shapiro.test(log10(meta$alpha))






## Heat map 1 #####
HM50 <-  as.data.frame(t(merged.pro.hell))
HM50 <- cbind(HM50, total = rowSums(HM50)) #calc total for each ASV
HM50 <- HM50[order(-HM50$total),] #order by most abundant
HM50 <- as.matrix(HM50[,1:971]) #remove total column
HM25 <- HM50[1:25,] #take only top 50 most abundant taxa
HM50 <- HM50[1:25,] #take only top 50 most abundant taxa

#prepare data for plot
rownames25 <- paste(fullmap.16s[rownames(HM25),7], " (", fullmap.16s[rownames(HM25),2],")",sep = "") #extract taxa names and proportion

rownames50 <- paste(fullmap.16s[rownames(HM50),7], " (", fullmap.16s[rownames(HM50),2],")",sep = "") #extract taxa names and proportion

#function to scale data
cal_z_score <- function(x){
  (x - min(x)) / (max(x)-min(x))
}

library(pheatmap)
# using the scaling colors 
HM50scaled <- t(apply(HM50, 1, cal_z_score))
HM25scaled <- t(apply(HM25, 1, cal_z_score))
pheatmap(HM25)
pheatmap(HM50)

#Italicized names 
fancynames25 <- lapply(
  rownames25,
  function(x) bquote(italic(.(x))))

library(vegan)
#Calulating Bray dist 
distrow25 <- vegdist(HM25, "bray", diag = TRUE)

library(RColorBrewer)
heat.col <- colorRampPalette(brewer.pal(9, "Blues"))(100)
ancols <- list("Taxonomic Mode" = cols2)
rownames(merged.pro.meta) <- merged.pro.meta$libraryname_pro

eep <- merged.pro.meta[colnames(HM25),]

colannotation <- data.frame("Taxonomic Mode" = eep$pro_mode_hell_k11)
#"Sample Type" = eep$type) #create annotation
rownames(colannotation) <- rownames(eep)

#25 most abundant sig ASVs
pdf("figs/mostabundanthm.pdf", width = 15, height = 10)
pheatmap(HM25scaled, labels_row = as.expression(fancynames25), annotation_col = colannotation, col = heat.col, show_colnames = F, cluster_rows = TRUE, cluster_cols = TRUE)
dev.off()


### Heat map 2 ####
library(pheatmap)


rownames(fullmap.16s)
colnames(merged.pro.hell) 

x <- rbind(map.arc,map.bac)


list.pro.modes <- c("promode1", "promode2", "promode3", "promode4", "promode5", "promode6", "promode7", "promode8", "promode9","promode10", "promode11")


map.pro.modes <- x[colnames(merged.pro.hell),]
map.pro.modes$mode1 <- NA
map.pro.modes$mode2 <- NA
map.pro.modes$mode3 <- NA
map.pro.modes$mode4 <- NA
map.pro.modes$mode5 <- NA
map.pro.modes$mode6 <- NA
map.pro.modes$mode7 <- NA
map.pro.modes$mode8 <- NA
map.pro.modes$mode9 <- NA
map.pro.modes$mode10 <- NA
map.pro.modes$mode11 <- NA

pro.top50.ASVs <- {}
pro.top30.ASVs <- {}
pro.top100.ASVs <- {}


for (i in 1:length(list.pro.modes)) {
  mode <- merged.pro.hell[rownames(merged.pro.meta[which(merged.pro.meta$pro_mode_hell_k11 == paste0(i)),]),]
  mode <- as.data.frame(t(mode))
  mode$Total <- rowSums(mode)
  map.pro.modes[,2+i] <- mode[rownames(map.pro.modes),ncol(mode)]/(ncol(mode)-1)
  mode <- mode[order(-mode$Total),]
  write.csv(cbind(mode,merged.pro.map[rownames(mode),]), paste0("~/mosaic/new_2023/all_SOM/SOM_composition/ProMode",i,".csv"))
  top3 <- rownames(mode)[1:3]
  top5 <- rownames(mode)[1:5]
  top10 <- rownames(mode)[1:10]
  pro.top30.ASVs <- c(pro.top30.ASVs, top3)
  pro.top50.ASVs <- c(pro.top50.ASVs, top5)
  pro.top100.ASVs <- c(pro.top100.ASVs, top10)
}


map.pro.modes$ASV <- rownames(map.pro.modes)
map.pro.modes$ASVdup <-gsub("\\..*","",map.pro.modes$ASV)  
map.pro.modesnodup <- map.pro.modes[!duplicated(map.pro.modes[,15]),]

HMtotals <- map.pro.modes[pro.top50.ASVs,]
HMtotalsnodup <- HMtotals[!duplicated(HMtotals[,15]),]

HMtotals <- HMtotalsnodup[,3:13]


#HMtotals <- map.pro.modes[pro.top50.ASVs,3:13]



#HMalltest <- t(merged.pro.RA[,pro.top50.ASVs])
#pheatmap(HMalltest, show_colnames = FALSE, show_rownames = FALSE)
pheatmap(HMtotals, show_colnames = FALSE, show_rownames = FALSE)
#going with totals

#prepare data for plot
rownames <- paste(merged.pro.map[rownames(HMtotals),7], " (", merged.pro.map[rownames(HMtotals),2],")",sep = "") #extract taxa names and proportion
rownames[c(8,23,30)] <- "Bacteria (1)"

#function to scale data
cal_z_score <- function(x){
  (x - min(x)) / (max(x)-min(x))
}


# using the scaling colors 
HMtotalsscale <- t(apply(HMtotals, 1, cal_z_score))


#Italicized names 
fancynames <- lapply(
  rownames,
  function(x) bquote(italic(.(x))))

library(RColorBrewer)
heat.col <- colorRampPalette(brewer.pal(9, "Blues"))(100)
ancols <- list("Taxonomic Mode" = cols2)

HM16S <- pheatmap(HMtotalsscale, labels_row = rownames, col = heat.col, show_colnames = T, cluster_rows = TRUE, cluster_cols = TRUE)

pdf("figs/16smodeHM.pdf", width = 10, height = 8)
HM16S
dev.off()

#### Heatmap pathways #####

pathways <- read.csv("~/mosaic/mosaic_16S_18S/data/round2/20220919_mosaic.bacteria.path_tally.csv", row.names = 1)
pathways <- pathways[rownames(merged.pro.hell),]
pathways <- as.data.frame(t(pathways))


pathways <- cbind(pathways, total = rowSums(pathways)) #calc total for each pathway
pathways <- pathways[order(-pathways$total),] #order by most abundant
pathways <- as.matrix(pathways[,1:971]) #remove total column
HM30 <- pathways[1:30,] #take only top 30 most abundant
HM15 <- pathways[1:15,] 



#function to scale data
cal_z_score <- function(x){
  (x - min(x)) / (max(x)-min(x))
}


# using the scaling colors 
HM15scaled <- t(apply(HM15, 1, cal_z_score))
HM30scaled <- t(apply(HM30, 1, cal_z_score))


distrow25 <- vegdist(HM30, "bray", diag = TRUE)

library(RColorBrewer)
heat.col <- colorRampPalette(brewer.pal(9, "Blues"))(100)
ancols <- list("Taxonomic Mode" = cols2)

eep <- merged.pro.meta[colnames(HM30),]

colannotation <- data.frame("Taxonomic Mode" = eep$pro_mode_hell_k11)
#"Sample Type" = eep$type) #create annotation
rownames(colannotation) <- rownames(eep)

#25 most abundant sig ASVs
pdf("figs/mostabundantpathways.pdf", width = 15, height = 10)
pheatmap(HM15, annotation_col = colannotation, col = heat.col, show_colnames = F, cluster_rows = TRUE, cluster_cols = TRUE)
dev.off()

### diveristy ######
bkleh <- as.data.frame(t(merged.pro.RA))
merged.pro.meta$alpha <- diversity(bkleh[,rownames(merged.pro.meta)],
                        MARGIN = 2,
                        index = "invsimpson")

shapiro.test(log10(merged.pro.meta$alpha))

kt <-  merged.pro.meta %>% kruskal_test(alpha ~ pro_mode_hell_k11)
pdf("figs/boxplots/diversity.pdf", width = 10, height = 5)
ggplot(merged.pro.meta, aes(x= pro_mode_hell_k11, y = alpha, color = pro_mode_hell_k11)) + geom_boxplot() + theme_light() + xlab("Taxonomic Mode") + ylab("Inverse Simpson Diversity Index") + theme(legend.position = "none",  legend.title = element_blank(), axis.title = element_text(size=16), axis.text = element_text(size = 14)) + scale_color_manual(values = cols2) + labs(subtitle = get_test_label(kt, detailed = TRUE))
dev.off()




###### create odv file ######



odvfile <- merged.pro.meta[,1:2]
odvfile$Cruise <- merged.pro.meta$type
odvfile$Station <- merged.pro.meta$dship_operation
odvfile$Date <- merged.pro.meta$date
odvfile$latitudedeg <- merged.pro.meta$lat
odvfile$longitudedeg <- merged.pro.meta$lon
odvfile$BotDept <- -999
odvfile$Niskin <- merged.pro.meta$niskin.bottle
odvfile$Depth <- merged.pro.meta$depth_m
odvfile$library <- merged.pro.meta$libraryname_pro

write.csv(odvfile, "odvplotting.csv")

