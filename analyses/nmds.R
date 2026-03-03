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
# Continuous traits
d$logSeedWeight <- log(d$seedWeights)
conifer$logSeedWeight <- log(conifer$seedWeights)
angio$logSeedWeight <- log(angio$seedWeights)

d$logFruit <- log(d$fruitSizeAve)
conifer$logFruit <- log(conifer$fruitSizeAve)
angio$logFruit <- log(angio$fruitSizeAve)

d$logSeedSize <- log(d$seedSizeAve)
conifer$logSeedSize <- log(conifer$seedSizeAve)
angio$logSeedSize <- log(angio$seedSizeAve)

phytree$node.label <- NULL

phyconifer <- drop.tip(phytree, setdiff(phytree$tip.label, conifer$latbi))
phyangio   <- drop.tip(phytree, setdiff(phytree$tip.label, angio$latbi))

## Make species x traits matrix ----
matrixCon <- data.frame(conifer$droughtTolerance,conifer$typeMonoOrDio,conifer$seedDispersal,conifer$seedDormancy,conifer$mastEvent,conifer$leafLongevity,conifer$oilContent,conifer$logSeedWeight,conifer$logFruit,conifer$logSeedSize,conifer$latbi)
colnames(matrixCon) <- c('Drought Tolerance','Reproductive Type','Dispersal Mode','Dormancy','Mast','Leaf Longevity','Oil Content','Seed Weight','Fruit Size','Seed Size','Species')
matrixConNoM <- data.frame(conifer$droughtTolerance,conifer$typeMonoOrDio,conifer$seedDispersal,conifer$seedDormancy,conifer$leafLongevity,conifer$oilContent,conifer$logSeedWeight,conifer$logFruit,conifer$logSeedSize,conifer$latbi)
colnames(matrixConNoM) <- c('Drought Tolerance','Reproductive Type','Dispersal Mode','Dormancy','Leaf Longevity','Oil Content','Seed Weight','Fruit Size','Seed Size','Species')


matrixAngio <- data.frame(angio$droughtTolerance,angio$typeMonoOrDio,angio$pollination,angio$seedDispersal,angio$seedDormancy,angio$mastEvent,angio$leafLongevity,angio$oilContent,angio$logSeedWeight,angio$logFruit,angio$logSeedSize,angio$latbi)
colnames(matrixAngio) <- c('Drought Tolerance','Reproductive Type','Pollination','Dispersal Mode','Dormancy','Mast','Leaf Longevity','Oil Content','Seed Weight','Fruit Size','Seed Size','Species')
#levels(matrixAngio$Pollination)[levels(matrixAngio$Pollination) == "wind and animals"] <- "both"
#matrixAngio$Pollination <- addNA(as.factor(matrixAngio$Pollination))
#levels(matrixAngio$Pollination)[is.na(levels(matrixAngio$Pollination))] <- "Missing"
matrixAngioNoM <- data.frame(angio$droughtTolerance,angio$typeMonoOrDio,angio$pollination,angio$seedDispersal,angio$seedDormancy, angio$leafLongevity,angio$oilContent,angio$logSeedWeight,angio$logFruit,angio$logSeedSize,angio$latbi)
colnames(matrixAngioNoM) <- c('Drought Tolerance','Reproductive Type','Pollination','Dispersal Mode','Dormancy','Leaf Longevity','Oil Content','Seed Weight','Fruit Size','Seed Size','Species')

rownames(matrixCon) <- matrixCon$Species
matrixCon$Species <- NULL

rownames(matrixAngio)   <- matrixAngio$Species
matrixAngio$Species <- NULL

rownames(matrixConNoM) <- matrixConNoM$Species
matrixConNoM$Species <- NULL

rownames(matrixAngioNoM)   <- matrixAngioNoM$Species
matrixAngioNoM$Species <- NULL
### Calculate the gower distance ----
#con_gower <- daisy(matrixCon, metric = "gower")
#angio_gower <- daisy(matrixAngio, metric = "gower")
con_gower_noM <- daisy(matrixConNoM, metric = "gower")
angio_gower_noM <- daisy(matrixAngioNoM, metric = "gower")

### NMDS ----
set.seed(11223)

nmds_con_noM <- metaMDS(
  con_gower_noM,
  k = 3,
  trymax = 200,
  autotransform = FALSE
)


conScoresN <- as.data.frame(scores(nmds_con_noM))
conScoresN$"Drought" <- matrixConNoM$"Drought Tolerance"
conScoresN$"Reproductive" <- matrixConNoM$"Reproductive Type"
conScoresN$"Dispersal" <- matrixConNoM$"Dispersal Mode"
conScoresN$"Dormancy" <- matrixConNoM$"Dormancy"
conScoresN$"Mast" <- matrixCon$"Mast"
conScoresN$"Weight" <- matrixCon$"Seed Weight"


set.seed(11223)


nmds_angio_N <- metaMDS(
  angio_gower_noM,
  k = 3,
  trymax = 100
)

#stressplot(nmds_angio)


angioScoresN <- as.data.frame(scores(nmds_angio_N))
angioScoresN$"Drought" <- matrixAngio$"Drought Tolerance"
angioScoresN$"Reproductive" <- matrixAngio$"Reproductive Type"
angioScoresN$"Dispersal" <- matrixAngio$"Dispersal Mode"
angioScoresN$"Dormancy" <- matrixAngio$"Dormancy"
angioScoresN$"Mast" <- matrixAngio$"Mast"
angioScoresN$"Weight" <- matrixAngio$"Seed Weight"


conN12 <- ggplot(conScoresN, aes(NMDS1, NMDS2)) +
  geom_point(aes(
    colour = Drought,
    fill   = interaction(Drought, Mast),
    shape  = Reproductive,
    size   = Weight,
    alpha  = factor(Mast)
  ),
  stroke = 1.2
  ) +
  scale_shape_manual(values = c(21, 22), guide = "none") +   # fillable shapes only
  scale_colour_manual(values = c("#95B958","#F4D166","#6194BF"),guide = "none") +
  scale_fill_manual(
    values = c(
      "Low.0"    = "white",
      "Moderate.0" = "white",
      "High.0"   = "white",
      "Low.1"   = "#F4D166",
      "Moderate.1"= "#6194BF",
      "High.1"  = "#95B958"
    ), guide = "none"
  ) + scale_alpha_manual(
    name = "Mast",
    values = c("0" = 1, "1" = 1),   # keep plot fully opaque
    guide = guide_legend(
      override.aes = list(
        shape = 21,
        size = 5,
        colour = "black",
        fill = c("white", "black")
      )
    )
  ) +
  scale_size_continuous(range = c(1, 10), guide = "none") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30"))



conN23 <- ggplot(conScoresN, aes(NMDS2, NMDS3)) +
  geom_point(aes(
    colour = Drought,
    fill   = interaction(Drought, Mast),
    shape  = Reproductive,
    size   = Weight,
    alpha  = factor(Mast)
  ),
  stroke = 1.2
  ) +
  scale_shape_manual(values = c(21, 22), guide = "none") +   # fillable shapes only
  scale_colour_manual(values = c("#95B958","#F4D166","#6194BF"),guide = "none") +
  scale_fill_manual(
    values = c(
      "Low.0"    = "white",
      "Moderate.0" = "white",
      "High.0"   = "white",
      "Low.1"   = "#F4D166",
      "Moderate.1"= "#6194BF",
      "High.1"  = "#95B958"
    ), guide = "none"
  ) + scale_alpha_manual(
    name = "Mast",
    values = c("0" = 1, "1" = 1),   # keep plot fully opaque
    guide = guide_legend(
      override.aes = list(
        shape = 21,
        size = 5,
        colour = "black",
        fill = c("white", "black")
      )
    )
  ) +
  scale_size_continuous(range = c(1, 10), guide = "none") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30"))

# Angiosperm masting yes and masting no


angioN12 <- ggplot(angioScoresN, aes(NMDS1, NMDS2)) +
  geom_point(aes(
    colour = Drought,
    fill   = interaction(Drought, Mast),
    shape  = Reproductive,
    size   = Weight,
    alpha  = factor(Mast)
  ),
  stroke = 1.2
  ) +
  scale_shape_manual(values = c(21, 22, 23), na.value = 24) +   # fillable shapes only
  scale_colour_manual(values = c("#95B958","#F4D166","#6194BF")) +
  scale_fill_manual(
    values = c(
      "Low.0"    = "white",
      "Moderate.0" = "white",
      "High.0"   = "white",
      "Low.1"   = "#F4D166",
      "Moderate.1"= "#6194BF",
      "High.1"  = "#95B958"
    )
  ) + scale_alpha_manual(
    name = "Mast",
    values = c("0" = 1, "1" = 1),   # keep plot fully opaque
    guide = guide_legend(
      override.aes = list(
        shape = 21,
        size = 5,
        colour = "black",
        fill = c("white", "black")
      )
    )
  ) +
  scale_size_continuous(range = c(1, 10)) +
  guides(fill = "none") + 
  labs(
    colour = "Drought Tolerance",
    shape  = "Reproductive Type",
    size   = "Seed Weight",
    alpha  = "Mast"
  ) +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30"))



angioN23 <- ggplot(angioScoresN, aes(NMDS2, NMDS3)) +
  geom_point(aes(
    colour = Drought,
    fill   = interaction(Drought, Mast),
    shape  = Reproductive,
    size   = Weight,
    alpha  = factor(Mast)
  ),
  stroke = 1.2
  ) +
  scale_shape_manual(values = c(21, 22, 23), na.value = 24) +   # fillable shapes only
  scale_colour_manual(values = c("#95B958","#F4D166","#6194BF")) +
  scale_fill_manual(
    values = c(
      "Low.0"    = "white",
      "Moderate.0" = "white",
      "High.0"   = "white",
      "Low.1"   = "#F4D166",
      "Moderate.1"= "#6194BF",
      "High.1"  = "#95B958"
    )
  ) + scale_alpha_manual(
    name = "Mast",
    values = c("0" = 1, "1" = 1),   # keep plot fully opaque
    guide = guide_legend(
      override.aes = list(
        shape = 21,
        size = 5,
        colour = "black",
        fill = c("white", "black")
      )
    )
  ) +
  scale_size_continuous(range = c(1, 10)) +
  guides(fill = "none") + 
  labs(
    colour = "Drought Tolerance",
    shape  = "Reproductive Type",
    size   = "Seed Weight",
    alpha  = "Mast"
  )  +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30"))

angio <- angioN12 + angioN23
conifer <- conN12 + conN23

final <- angio / conifer +
  plot_layout(heights = c(1, 1), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    legend.position = "right",
    plot.tag = element_text(face = "bold", size = 8)
  )

ggsave(
  filename = "output/figures/nmdsFinal.pdf",
  plot = final,
  width = 12,
  height = 10
)

phylomorphospace(
  phyconifer,
  nmds_con$points,
  label = "horizontal"
)





