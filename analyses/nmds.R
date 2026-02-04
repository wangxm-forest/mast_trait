#### Phylogenetic PCA ####
## Started by Mao ##
## Jan-16-2026 ##

library(ape)
library(phytools)
library(cluster)
library(vegan)
library(plotly)
library(gridExtra)


rm(list = ls())
options(stringsAsFactors = FALSE)

setwd("C:/PhD/Project/PhD_thesis/mast_trait")

# Extract legend from p_all
get_legend <- function(myplot) {
  g <- ggplotGrob(myplot)
  leg <- g$grobs[which(sapply(g$grobs, function(x) x$name) == "guide-box")]
  leg[[1]]
}

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

matrixAngio <- data.frame(angio$droughtTolerance,angio$typeMonoOrDio,angio$pollination,angio$seedDispersal,angio$seedDormancy,angio$mastEvent,angio$leafLongevity,angio$oilContent,angio$logSeedWeightStd,angio$logFruitStd,angio$logSeedSizeStd,angio$latbi)
colnames(matrixAngio) <- c('Drought Tolerance','Reproductive Type','Pollination','Dispersal Mode','Dormancy','Mast','Leaf Longevity','Oil Content','Seed Weight','Fruit Size','Seed Size','Species')
#levels(matrixAngio$Pollination)[levels(matrixAngio$Pollination) == "wind and animals"] <- "both"
#matrixAngio$Pollination <- addNA(as.factor(matrixAngio$Pollination))
#levels(matrixAngio$Pollination)[is.na(levels(matrixAngio$Pollination))] <- "Missing"

AngioNum <- c("Leaf Longevity","Oil Content","Seed Weight","Fruit Size","Seed Size")
#AngioCat <- c("Drought Tolerance","Pollination","Reproductive Type","Dispersal Mode","Dormancy","Mast")
AngioNumD <- matrixAngio[AngioNum]
#AngioCatD <- matrixAngio[AngioCat]

rownames(matrixCon) <- matrixCon$Species
matrixCon$Species <- NULL

rownames(matrixAngio)   <- matrixAngio$Species
matrixAngio$Species <- NULL
### Calculate the gower distance ----
con_gower <- daisy(matrixCon, metric = "gower")
angio_gower <- daisy(matrixAngio, metric = "gower")

### NMDS ----
set.seed(11223)
nmds_con <- metaMDS(
  con_gower,
  k = 3,
  trymax = 200,
  autotransform = FALSE
)


nmds_con$stress
#nmds_con_vectors <- envfit(nmds_con, matrixCon, choices = c(1,2,3),permutations = 999, na.rm = TRUE)


conScores <- as.data.frame(scores(nmds_con))
conScores$"Drought" <- matrixCon$"Drought Tolerance"
conScores$"Reproductive" <- matrixCon$"Reproductive Type"
conScores$"Dispersal" <- matrixCon$"Dispersal Mode"
conScores$"Dormancy" <- matrixCon$"Dormancy"
conScores$"Mast" <- matrixCon$"Mast"


#enCoordCon <- as.data.frame(scores(nmds_con_vectors, "vectors")) * ordiArrowMul(nmds_con_vectors)

# Multiplier for the first plot
#enCoordCon12 <- as.data.frame(scores(nmds_con_vectors, "vectors", choices = c(1, 2))) * ordiArrowMul(nmds_con_vectors, choices = c(1, 2))
# Multiplier for the second plot
#enCoordCon23 <- as.data.frame(scores(nmds_con_vectors, "vectors", choices = c(2, 3))) * ordiArrowMul(nmds_con_vectors, choices = c(2, 3))


#plot_ly(conScores, x = ~NMDS1, y = ~NMDS2, z = ~NMDS3, 
#        color = ~Mast, type = 'scatter3d', mode = 'markers')

conDrought12 <- ggplot(data = conScores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = Drought, shape = Mast), size = 3, alpha = 0.6) + scale_shape_manual(values = c(1, 16)) +
  scale_colour_manual(values = c("orange","steelblue", "darkgreen"), 
                      na.value = "grey80") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Drought Tolerance", shape = "Mast")

conDrought23 <- ggplot(data = conScores, aes(x = NMDS2, y = NMDS3)) + 
  geom_point(aes(colour = Drought, shape = Mast), size = 3, alpha = 0.6) + scale_shape_manual(values = c(1, 16)) +
  scale_colour_manual(values = c("orange","steelblue", "darkgreen"), 
                      na.value = "grey80") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Drought Tolerance", shape = "Mast")

shared_legend1 <- get_legend(conDrought12)

p1 <- conDrought12 + theme(legend.position = "none")
p2 <- conDrought23 + theme(legend.position = "none")


conRep12 <- ggplot(data = conScores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = Reproductive, shape = Mast), size = 3, alpha = 0.6) + scale_shape_manual(values = c(1, 16)) +
  scale_colour_manual(values = c("orange","steelblue", "darkgreen"), 
                      na.value = "grey80") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Reproductive Type", shape = "Mast")



conRep23 <- ggplot(data = conScores, aes(x = NMDS2, y = NMDS3)) + 
  geom_point(aes(colour = Reproductive, shape = Mast), size = 3, alpha = 0.6) + scale_shape_manual(values = c(1, 16)) +
  scale_colour_manual(values = c("orange","steelblue", "darkgreen"), 
                      na.value = "grey80") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Reproductive Type", shape = "Mast")

shared_legend2 <- get_legend(conRep12)
p3 <- conRep12 + theme(legend.position = "none")
p4 <- conRep23 + theme(legend.position = "none")

conDis12 <- ggplot(data = conScores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = Dispersal, shape = Mast), size = 3, alpha = 0.6) + scale_shape_manual(values = c(1, 16)) +
  scale_colour_manual(values = c("orange","steelblue", "darkgreen"), 
                      na.value = "grey80") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Dispersal Type", shape = "Mast")

conDis23 <- ggplot(data = conScores, aes(x = NMDS2, y = NMDS3)) + 
  geom_point(aes(colour = Dispersal, shape = Mast), size = 3, alpha = 0.6) + scale_shape_manual(values = c(1, 16)) +
  scale_colour_manual(values = c("orange","steelblue", "darkgreen"), 
                      na.value = "grey80") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Dispersal Type", shape = "Mast")

shared_legend3 <- get_legend(conDis12)
p5 <- conDis12 + theme(legend.position = "none")
p6 <- conDis23 + theme(legend.position = "none")


conDor12 <- ggplot(data = conScores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = Dormancy, shape = Mast), size = 3, alpha = 0.6) + scale_shape_manual(values = c(1, 16)) +
  scale_colour_manual(values = c("orange","steelblue"), 
                      na.value = "grey80") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Dormancy", shape = "Mast")

conDor23 <- ggplot(data = conScores, aes(x = NMDS2, y = NMDS3)) + 
  geom_point(aes(colour = Dormancy, shape = Mast), size = 3, alpha = 0.6) + scale_shape_manual(values = c(1, 16)) +
  scale_colour_manual(values = c("orange","steelblue"), 
                      na.value = "grey80") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Dormancy", shape = "Mast")

shared_legend4 <- get_legend(conDor12)
p7 <- conDor12 + theme(legend.position = "none")
p8 <- conDor23 + theme(legend.position = "none")

pdf("output/figures/nmdsCon.pdf", width = 13, height = 20)

grid.arrange(p1, p2, shared_legend1, p3, p4, shared_legend2, p5, p6, shared_legend3, p7, p8, shared_legend4, 
             ncol = 3, 
             widths = c(3, 3, 0.8))

dev.off()


set.seed(11223)
nmds_angio <- metaMDS(
  angio_gower,
  k = 3,
  trymax = 100
)

#stressplot(nmds_angio)
nmds_angio$stress

angioScores <- as.data.frame(scores(nmds_angio))

angioScores$"Drought" <- matrixAngio$"Drought Tolerance"
angioScores$"Reproductive" <- matrixAngio$"Reproductive Type"
angioScores$"Dispersal" <- matrixAngio$"Dispersal Mode"
angioScores$"Dormancy" <- matrixAngio$"Dormancy"
angioScores$"Pollination" <- matrixAngio$"Pollination"
angioScores$"Mast" <- matrixAngio$"Mast"

angioDrought12 <- ggplot(data = angioScores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = Drought, shape = Mast), size = 3, alpha = 0.6) + scale_shape_manual(values = c(1, 16)) +
  scale_colour_manual(values = c("orange","steelblue", "darkgreen"), 
                      na.value = "grey80") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Drought Tolerance", shape = "Mast")

angioDrought23 <- ggplot(data = angioScores, aes(x = NMDS2, y = NMDS3)) + 
  geom_point(aes(colour = Drought, shape = Mast), size = 3, alpha = 0.6) + scale_shape_manual(values = c(1, 16)) +
  scale_colour_manual(values = c("orange","steelblue", "darkgreen"), 
                      na.value = "grey80") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Drought Tolerance", shape = "Mast")

shared_legend1 <- get_legend(angioDrought12)
p1 <- angioDrought12 + theme(legend.position = "none")
p2 <- angioDrought12 + theme(legend.position = "none")

angioRep12 <- ggplot(data = angioScores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = Reproductive, shape = Mast), size = 3, alpha = 0.6) + scale_shape_manual(values = c(1, 16)) +
  scale_colour_manual(values = c("orange","steelblue", "darkgreen"), 
                      na.value = "grey80") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Reproductive Type", shape = "Mast")

angioRep23 <- ggplot(data = angioScores, aes(x = NMDS2, y = NMDS3)) + 
  geom_point(aes(colour = Reproductive, shape = Mast), size = 3, alpha = 0.6) + scale_shape_manual(values = c(1, 16)) +
  scale_colour_manual(values = c("orange","steelblue", "darkgreen"), 
                      na.value = "grey80") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Reproductive Type", shape = "Mast")

shared_legend2 <- get_legend(angioRep12)
p3 <- angioRep12 + theme(legend.position = "none")
p4 <- angioRep23 + theme(legend.position = "none")

angioDis12 <- ggplot(data = angioScores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = Dispersal, shape = Mast), size = 3, alpha = 0.6) + scale_shape_manual(values = c(1, 16)) +
  scale_colour_manual(values = c("orange","steelblue", "darkgreen"), 
                      na.value = "grey80") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Dispersal Type", shape = "Mast")

angioDis23 <- ggplot(data = angioScores, aes(x = NMDS2, y = NMDS3)) + 
  geom_point(aes(colour = Dispersal, shape = Mast), size = 3, alpha = 0.6) + scale_shape_manual(values = c(1, 16)) +
  scale_colour_manual(values = c("orange","steelblue", "darkgreen"), 
                      na.value = "grey80") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Dispersal Type", shape = "Mast")

shared_legend3 <- get_legend(angioDis12)
p5 <- angioDis12 + theme(legend.position = "none")
p6 <- angioDis23 + theme(legend.position = "none")

angioDor12 <- ggplot(data = angioScores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = Dormancy, shape = Mast), size = 3, alpha = 0.6) + scale_shape_manual(values = c(1, 16)) +
  scale_colour_manual(values = c("orange","steelblue"), 
                      na.value = "grey80") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Dormancy", shape = "Mast")

angioDor23 <- ggplot(data = angioScores, aes(x = NMDS2, y = NMDS3)) + 
  geom_point(aes(colour = Dormancy, shape = Mast), size = 3, alpha = 0.6) + scale_shape_manual(values = c(1, 16)) +
  scale_colour_manual(values = c("orange","steelblue"), 
                      na.value = "grey80") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Dormancy", shape = "Mast")


shared_legend4 <- get_legend(angioDor12)
p7 <- angioDor12 + theme(legend.position = "none")
p8 <- angioDor23 + theme(legend.position = "none")

angioPoll12 <- ggplot(data = angioScores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = Pollination, shape = Mast), size = 3, alpha = 0.6) + scale_shape_manual(values = c(1, 16)) +
  scale_colour_manual(values = c("orange","steelblue", "darkgreen"), 
                      na.value = "grey80") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Pollination", shape = "Mast")

angioPoll23 <- ggplot(data = angioScores, aes(x = NMDS2, y = NMDS3)) + 
  geom_point(aes(colour = Pollination, shape = Mast), size = 3, alpha = 0.6) + scale_shape_manual(values = c(1, 16)) +
  scale_colour_manual(values = c("orange","steelblue", "darkgreen"), 
                      na.value = "grey80") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Pollination", shape = "Mast")

shared_legend5 <- get_legend(angioPoll12)
p9 <- angioPoll12 + theme(legend.position = "none")
p10 <- angioPoll23 + theme(legend.position = "none")

pdf("output/figures/nmdsAngio.pdf", width = 13, height = 25)

grid.arrange(p1, p2, shared_legend1, p3, p4, shared_legend2, p5, p6, shared_legend3, p7, p8, shared_legend4, p9, p10, shared_legend5,
             ncol = 3, 
             widths = c(3, 3, 0.8))

dev.off()

phylomorphospace(
  phyconifer,
  nmds_con$points,
  label = "horizontal"
)





