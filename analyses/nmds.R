#### Phylogenetic PCA ####
## Started by Mao ##
## Jan-16-2026 ##

library(ape)
library(phytools)
library(cluster)
library(vegan)
library(plotly)


rm(list = ls())
options(stringsAsFactors = FALSE)

setwd("C:/PhD/Project/PhD_thesis/mast_trait")

# Load & prepare Data
d <- read.csv("data/cleanSilvics.csv")
phytree <- read.tree("output/silvicsPhylogenyFull.tre")

d$latbi <- gsub(" ", "_", d$latbi)
d <- d[!is.na(d$mastEvent), ]
d$mastEvent <- ifelse(d$mastEvent == "Y", 1, 0)

# Create group variable
d$group <- ifelse(d$familyName %in% c("Pinaceae", "Taxodiaceae"),
                  "conifer", "angiosperm")
d$group <- factor(d$group)

conifer <- d[d$familyName %in% c("Pinaceae", "Taxodiaceae"), ]
angio   <- d[!(d$familyName %in% c("Pinaceae", "Taxodiaceae")), ]


# Factorize all the categorical traits, I am doing this separately because conifer and angio have different levels for some traits
conifer$droughtTolerance <- as.factor(conifer$droughtTolerance)
conifer$typeMonoOrDio <- as.factor(conifer$typeMonoOrDio)
conifer$pollination <- as.factor(conifer$pollination)
conifer$seedDispersal <- as.factor(conifer$seedDispersal)
conifer$seedDormancy <- as.factor(conifer$seedDormancy)
conifer$mastEvent <- as.factor(conifer$mastEvent)
angio$droughtTolerance <- as.factor(angio$droughtTolerance)
angio$typeMonoOrDio <- as.factor(angio$typeMonoOrDio)
angio$pollination <- as.factor(angio$pollination)
angio$seedDispersal <- as.factor(angio$seedDispersal)
angio$seedDormancy <- as.factor(angio$seedDormancy)
angio$mastEvent <- as.factor(angio$mastEvent)
# log10-transform + scale
log_scale <- function(x) scale(log10(x))

# Continuous traits
d$logSeedWeightStd <- log_scale(d$seedWeights)
conifer$logSeedWeightStd <- log_scale(conifer$seedWeights)
angio$logSeedWeightStd   <- log_scale(angio$seedWeights)

d$logFruitStd <- log_scale(d$fruitSizeAve)
conifer$logFruitStd <- log_scale(conifer$fruitSizeAve)
angio$logFruitStd   <- log_scale(angio$fruitSizeAve)

d$logSeedSizeStd <- log_scale(d$seedSizeAve)
conifer$logSeedSizeStd <- log_scale(conifer$seedSizeAve)
angio$logSeedSizeStd   <- log_scale(angio$seedSizeAve)

phytree$node.label <- NULL

phyconifer <- drop.tip(phytree, setdiff(phytree$tip.label, conifer$latbi))
phyangio   <- drop.tip(phytree, setdiff(phytree$tip.label, angio$latbi))

## Make species x traits matrix ----
matrixCon <- data.frame(conifer$droughtTolerance,conifer$typeMonoOrDio,conifer$seedDispersal,conifer$seedDormancy,conifer$mastEvent,conifer$leafLongevity,conifer$oilContent,conifer$logSeedWeightStd,conifer$logFruitStd,conifer$logSeedSizeStd,conifer$latbi)
colnames(matrixCon) <- c('Drought Tolerance','Reproductive Type','Dispersal Mode','Dormancy','Mast','Leaf Longevity','Oil Content','Seed Weight','Fruit Size','Seed Size','Species')

ConNum <- c("Leaf Longevity","Oil Content","Seed Weight","Fruit Size","Seed Size")
ConCat <- c("Drought Tolerance","Reproductive Type","Dispersal Mode","Dormancy","Mast")
ConNumD <- matrixCon[ConNum]
ConCatD <- matrixCon[ConCat]

matrixAngio <- data.frame(angio$droughtTolerance,angio$typeMonoOrDio,angio$seedDispersal,angio$seedDormancy,angio$mastEvent,angio$leafLongevity,angio$oilContent,angio$logSeedWeightStd,angio$logFruitStd,angio$logSeedSizeStd,angio$latbi)
colnames(matrixAngio) <- c('Drought Tolerance','Reproductive Type','Dispersal Mode','Dormancy','Mast','Leaf Longevity','Oil Content','Seed Weight','Fruit Size','Seed Size','Species')
#levels(matrixAngio$Pollination)[levels(matrixAngio$Pollination) == "wind and animals"] <- "both"
#matrixAngio$Pollination <- addNA(as.factor(matrixAngio$Pollination))
#levels(matrixAngio$Pollination)[is.na(levels(matrixAngio$Pollination))] <- "Missing"

#AngioNum <- c("Leaf Longevity","Oil Content","Seed Weight","Fruit Size","Seed Size")
#AngioCat <- c("Drought Tolerance","Pollination","Reproductive Type","Dispersal Mode","Dormancy","Mast")
#AngioNumD <- matrixAngio[AngioNum]
#AngioCatD <- matrixAngio[AngioCat]

rownames(matrixCon) <- matrixCon$Species
matrixCon$Species <- NULL

rownames(matrixAngio)   <- matrixAngio$Species
matrixAngio$Species <- NULL
### Calculate the gower distance ----
con_gower <- daisy(matrixCon, metric = "gower")
angio_gower <- daisy(matrixAngio, metric = "gower")

### NMDS ----

nmds_con <- metaMDS(
  con_gower,
  k = 3,
  trymax = 200
)


nmds_con$stress
nmds_con_vectors <- envfit(nmds_con, matrixCon, permutations = 999, na.rm = TRUE)
nmds_con_vectors



plot(nmds_con, type = "n") 
points(nmds_con, pch = 19)
plot(nmds_con_vectors, cex = 0.8)

conScores <- as.data.frame(scores(nmds_con))
conScores$"Drought" <- matrixCon$"Drought Tolerance"
conScores$"Reproductive" <- matrixCon$"Reproductive Type"
conScores$"Dispersal" <- matrixCon$"Dispersal Mode"
conScores$"Dormancy" <- matrixCon$"Dormancy"
conScores$"Mast" <- matrixCon$"Mast"


enCoordCon <- as.data.frame(scores(nmds_con_vectors, "vectors")) * ordiArrowMul(nmds_con_vectors)

plot_ly(conScores, x = ~NMDS1, y = ~NMDS2, z = ~NMDS3, 
        color = ~Mast, type = 'scatter3d', mode = 'markers')

conDrought1 <- ggplot() + 
  geom_point(data = conScores, aes(x = NMDS1, y = NMDS2, colour = Mast, shape = Drought), size = 3, alpha = 0.5) + scale_shape_manual(values = c(16, 17, 18, 19)) +
  scale_colour_manual(values = c("orange", "steelblue"))  + 
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = enCoordCon, linewidth = 1, alpha = 0.5, colour = "grey30") +
  geom_text(data = enCoordCon, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(enCoordCon))  +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Mast", shape = "Drought Tolerance")

conDrought2 <- ggplot() + 
  geom_point(data = conScores, aes(x = NMDS1, y = NMDS3, colour = Mast, shape = Drought), size = 3, alpha = 0.5) + scale_shape_manual(values = c(16, 17, 18, 19)) +
  scale_colour_manual(values = c("orange", "steelblue"))  + 
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS3), 
               data = enCoordCon, linewidth = 1, alpha = 0.5, colour = "grey30") +
  geom_text(data = enCoordCon, aes(x = NMDS1, y = NMDS3), colour = "grey30", 
            fontface = "bold", label = row.names(enCoordCon))  +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Mast", shape = "Drought Tolerance")

conDrought1
conDrought2

conReproductive <- ggplot() + 
  geom_point(data = conScores, aes(x = NMDS1, y = NMDS2, colour = Mast, shape = Reproductive), size = 3, alpha = 0.5) + scale_shape_manual(values = c(16, 17, 18, 19)) +
  scale_colour_manual(values = c("orange", "steelblue"))  + 
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = enCoordCon, linewidth = 1, alpha = 0.5, colour = "grey30") +
  geom_text(data = enCoordCon, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(enCoordCon)) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Mast", shape = "Reproductive type")

conReproductive

nmds_angio <- metaMDS(
  angio_gower,
  k = 3,
  trymax = 200
)

nmds_angio$stress
nmds_angio_vectors <- envfit(nmds_angio, matrixAngio, permutations = 999, na.rm = TRUE)

all(rownames(nmds_angio) == rownames(matrixAngio))
conScores <- as.data.frame(scores(nmds_con))
nmds_angio_vectors
plot(nmds_angio, type = "n") 
points(nmds_angio, pch = 19)
plot(nmds_angio, cex = 0.8)

angioScores <- as.data.frame(scores(nmds_angio))
angioScores$"Mast" <- matrixAngio$"Mast"
plot_ly(angioScores, x = ~NMDS1, y = ~NMDS2, z = ~NMDS3, 
        color = ~Mast, type = 'scatter3d', mode = 'markers')

angioScores$"Drought" <- matrixAngio$"Drought Tolerance"
angioScores$"Reproductive" <- matrixAngio$"Reproductive Type"
angioScores$"Dispersal" <- matrixAngio$"Dispersal Mode"
angioScores$"Dormancy" <- matrixAngio$"Dormancy"
angioScores$"Mast" <- matrixAngio$"Mast"


enCoordAngio <- as.data.frame(scores(nmds_angio_vectors, "vectors")) * ordiArrowMul(nmds_angio_vectors)

angioDrought <- ggplot() + 
  geom_point(data = angioScores, aes(x = NMDS1, y = NMDS2, colour = Mast, shape = Drought), size = 3, alpha = 0.5) + scale_shape_manual(values = c(16, 17, 18, 19)) +
  scale_colour_manual(values = c("orange", "steelblue"))  + 
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = enCoordAngio, linewidth = 1, alpha = 0.5, colour = "grey30") +
  geom_text(data = enCoordAngio, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(enCoordAngio))  +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Mast", shape = "Drought Tolerance")

angioDrought

phylomorphospace(
  phyconifer,
  nmds_con$points,
  label = "horizontal"
)


scores <- nmds_con$points   # rows = species


