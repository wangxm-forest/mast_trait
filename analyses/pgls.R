####PGLS####
##started by Mao##
##Nov-28-2025##

library(caper)
library(geiger)

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

# working directory
setwd("C:/PhD/Project/PhD_thesis/mast_trait")

# read in the clean data and phylogeny tree
d <- read.csv("data/cleanSilvics.csv", header = TRUE)
phytree <- read.tree("output/silvicsPhylogenyFull.tre")
d$latbi <- gsub(" ", "_", d$latbi)
# remove rows with mastEvent data being NA
d <- d[!is.na(d$mastEvent), ]
# Make a subset of conifers only
conifer <- subset(d, familyName %in% c("Pinaceae","Taxodiaceae"))
# Make a subset of angiosperm only
angio <- subset(d, !(familyName %in% c("Pinaceae","Taxodiaceae")))

# Prune the tree for conifer species only
phyconifer<-drop.tip(phytree, setdiff(phytree$tip.label, angio$latbi))
phyangio <- drop.tip(phytree,phytree$tip.label[-match(conifer$latbi, phytree$tip.label)])

# names match between the tree and the data frame
name.check(conifer, phyconifer)
name.check(angio, phyangio)

# Plot the tree
plot(phyconifer,,cex=0.4)
########################Mast (Y/N)###########################################
####
pglsModel2 <- gls(hostility ~ ecomorph, correlation = corBrownian(phy = anoleTree),
                  data = anoleData, method = "ML")