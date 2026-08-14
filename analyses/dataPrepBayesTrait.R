## Started 13 Aug 2026 ##
## Started by Mao ##
## Explore the BayesTraits for my trait analysis

library(phytools)
library(ape)
library(bayesTraitR)

setwd("C:/PhD/Project/PhD_thesis/mast_trait/")
d <- read.csv("data/cleanSilvicsSp.csv")

silvicsTree <- read.tree("output/silvicsPhylogenyFull.tre")
masting <- d[!duplicated(d$latbi), c("latbi","mastEvent")]
masting$mastEvent[is.na(masting$mastEvent)] <- "No information"
d$latbi <- gsub(" ", "_", d$latbi)
conifer <- d[d$familyName %in% c("Pinaceae", "Taxodiaceae"), ]
angio   <- d[!(d$familyName %in% c("Pinaceae", "Taxodiaceae")), ]


phyconifer <- drop.tip(silvicsTree, setdiff(silvicsTree$tip.label, conifer$latbi))
phyangio   <- drop.tip(silvicsTree, setdiff(silvicsTree$tip.label, angio$latbi))
d$logSeedWeight <- log(d$seedWeights)
conifer$logSeedWeight <- log(conifer$seedWeights)
angio$logSeedWeight   <- log(angio$seedWeights)

d$logFruit <- log(d$fruitSizeAve)
conifer$logFruit <- log(conifer$fruitSizeAve)
angio$logFruit <- log(angio$fruitSizeAve)

d$logSeedSize <- log(d$seedSizeAve)
conifer$logSeedSize <- log(conifer$seedSizeAve)
angio$logSeedSize   <- log(angio$seedSizeAve)


allspp <- silvicsTree$tip.label

angio$mastEvent[angio$mastEvent == "Y"]  <- 1
angio$mastEvent[angio$mastEvent == "N"]  <- 0

colnames(angio)
angio_mast <- subset(angio, 
                      select = c(latbi, mastEvent))

# Remove the underscores from species name since we need to write out the dataframe in txt file
angio_mast <- angio_mast[complete.cases(angio_mast), ]
angio_mast$latbi <- gsub("_", "", angio_mast$latbi)
write.table(angio_mast,"input/angio_mast.txt",
            col.names=F, row.names=F, quote=F,sep=" ")

# Do the same for the tree species labels
phyangio$tip.label <- gsub("_", "", phyangio$tip.label)

# BayesTrait requires the tree to be nexus format
writeNexus(phyangio,"output/phyangioNoSpace.tre")


setdiff(angio_mast$latbi,phyangio$tip.label)
setdiff(silvicsSpliced$tip.label, d$latbi)
