#### Phylogeny####
####Start by Mao####
####Sep 17####

library(ape)
library(phytools)
library(viridisLite)
library(geiger)
library(xtable)
library(stringr)
library(WorldFlora)
library(caper)
# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

setwd("C:/PhD/Project/PhD_thesis/mast_trait/")
d <- read.csv("data/cleanSilvics.csv")

maketree <- FALSE
if(maketree){
phy.plants<-read.tree("C:/PhD/Project/egret/analyses/input/ALLMB.tre")

backbone <- read.csv("C:/PhD/Project/classification.csv",head = TRUE, sep="\t")
d$speciesName <- tolower(d$speciesName)
# Remove trailing spaces:
d$speciesName <- str_trim(d$speciesName)
d$genusName <- str_trim(d$genusName)
d_species <- unique(paste(d$genusName, d$speciesName))
checks<-WFO.match(spec.data=d_species, WFO.data=backbone, counter=1, verbose=TRUE)
d_species_fix <- unique(checks$scientificName)
names_changed <- setdiff(d_species, d_species_fix)
names_changed

# Fix species names that Worldflora flagged #
d$genusName[which(d$genusName == "Chamaecyparis" & d$speciesName == "nootkatensis")] <- "Callitropsis"
d$speciesName[which(d$genusName == "Juniperus" & d$speciesName == "silicicola")] <- "virginiana_var._silicicola"
d$genusName[which(d$genusName == "Libocedrus" & d$speciesName == "decurrens")] <- "Calocedrus"
d$speciesName[which(d$genusName == "Acer" & d$speciesName == "barbatum")] <- "glabrum_var._glabrum"
d$speciesName[which(d$genusName == "Acer" & d$speciesName == "nigrum")] <- "saccharum_subsp._nigrum"
d$speciesName[which(d$genusName == "Aesculus" & d$speciesName == "octandra")] <- "flava"
d$speciesName[which(d$genusName == "Carya" & d$speciesName == "illinoensis")] <- "illinoinensis"
d$genusName[which(d$genusName == "Castanopsis" & d$speciesName == "chrysophylla")] <- "Chrysolepis"
d$genusName[which(d$genusName == "Lithocarpus" & d$speciesName == "densiflorus")] <- "Notholithocarpus"
d$speciesName[which(d$genusName == "Nyssa" & d$speciesName == "sylvatica var. biflora")] <- "biflora"
d$genusName[which(d$genusName == "Pithecellobium" & d$speciesName == "saman")] <- "Samanea"
d$speciesName[which(d$genusName == "Populus" & d$speciesName == "deltoides var occidentalis")] <- "deltoides_subsp._monilifera"
d$genusName[which(d$genusName == "Prosopis" & d$speciesName == "pallida")] <- "Neltuma"
d$speciesName[which(d$genusName == "Quercus" & d$speciesName == "falcata var. falcata")] <- "falcata"
d$speciesName[which(d$genusName == "Quercus" & d$speciesName == "falcata var. pagodifolia")] <- "pagoda"
d$speciesName[which(d$genusName == "Quercus" & d$speciesName == "nuttallii")] <- "texana"
d$latbi <- paste(d$genusName, d$speciesName, sep = "_")
unique(d$latbi)
d$latbi <- gsub(" ", "_", d$latbi)
## getting a list of genera in S&B's phylo
phy.genera<-unlist(
  lapply(strsplit(phy.plants$tip.label, "_"),function(x){return(x[1])})
)
phy.genera.uniq<-sort(unique(phy.genera))
phy.sps.uniqu <- phy.plants$tip.label
# just a check to confirm there are all the same format 
phy.sps.uniqu <- gsub(" ", "_", phy.sps.uniqu)

silvics.sps <- sort(unique(d$latbi))

silvics.phenosp.sps.inphylo <- silvics.sps[which(!silvics.sps%in%phy.sps.uniqu)] #15 out of 190

kew <- read.csv("C:/PhD/Project/egret/wcvp_Mao/wcvp_names.csv", header = TRUE, sep = "|", stringsAsFactors = FALSE )
kew$taxon_name <- gsub(" ", "_", kew$taxon_name)
kewvec <- unique(kew$taxon_name)
matchestree <- silvics.phenosp.sps.inphylo[silvics.phenosp.sps.inphylo %in% kewvec] 
matchestree

sppsilvicskew <- subset(d, latbi %in% silvics.phenosp.sps.inphylo)

sppsilvicskew2 <- sppsilvicskew[, c("latbi", "genusName","speciesName")]
kewsub <- subset(kew, taxon_name %in% silvics.phenosp.sps.inphylo)
# grab all parent ID for these species
accparentIDs <- kewsub$accepted_plant_name_id
# pull a vector of species names that correspond to these accepted IDs.
sub <- subset(kew, accepted_plant_name_id %in% accparentIDs)
# remove NAs
suby <- sub[sub$accepted_plant_name_id != "", ]
# pull a vector of all of these species 
kewnames <- suby$taxon_name
# look how many of these species are in the phylogeny tree:
withaccepted<-kewnames[which(kewnames%in%phy.sps.uniqu)] 
# subset kew for these species
matchednames <- subset(kew, taxon_name %in% withaccepted)

# remove unecessary columns
matchednamessub <- matchednames[, c("accepted_plant_name_id", "taxon_name")]

# add parent name ID in the df containing ALL the species we don't have match for in the tree
latbiwithId <- merge(sppsilvicskew2, kewsub[, c("taxon_name", "accepted_plant_name_id")], 
                     by.x = "latbi", by.y = "taxon_name", 
                     all.x = TRUE)

# merge matchednamessub and sppegretkew by accepted_plant_name_id
matchednamessilvics <- merge(latbiwithId, matchednamessub, by = "accepted_plant_name_id", all.x = TRUE) 
# change colnames
colnames(matchednamessilvics) <- c("accepted_plant_name_id", "silvicsname", "genusName", "speciesName", "sppMatch")
# this is the names that can be matched by code to kew's data base. 
head(matchednamessilvics)
nomatch <- matchednamessilvics[which(is.na(matchednamessilvics$sppMatch)),]
nrow(nomatch)

#phy.sps.uniqu[grepl("Magnolia_acuminata", phy.sps.uniqu)]
matchednamessilvics$sppMatch[which(matchednamessilvics$silvicsname == "Magnolia_acuminata")] <-  "Magnolia_acuminata_var._acuminata"

#phy.sps.uniqu[grepl("Magnolia_fraseri", phy.sps.uniqu)]
matchednamessilvics$sppMatch[which(matchednamessilvics$silvicsname == "Magnolia_fraseri")] <-  "Magnolia_fraseri_var._fraseri"

#phy.sps.uniqu[grepl("Metrosideros_polymorpha", phy.sps.uniqu)]
matchednamessilvics$sppMatch[which(matchednamessilvics$silvicsname == "Metrosideros_polymorpha")] <-  "Metrosideros_polymorpha_var._glaberrima"

#phy.sps.uniqu[grepl("Quercus_texana", phy.sps.uniqu)]
matchednamessilvics$sppMatch[which(matchednamessilvics$silvicsname == "Quercus_texana")] <-  "Quercus_buckleyi"

#phy.sps.uniqu[grepl("Notholithocarpus_densiflorus", phy.sps.uniqu)]
matchednamessilvics$sppMatch[which(matchednamessilvics$silvicsname == "Notholithocarpus_densiflorus")] <-  "Notholithocarpus_densiflorus_var._densiflorus"

#phy.sps.uniqu[grepl("Calophyllum_calaba", phy.sps.uniqu)]
matchednamessilvics$sppMatch[which(matchednamessilvics$silvicsname == "Calophyllum_calaba")] <-  "Calophyllum_calaba_var._bracteatum"

#phy.sps.uniqu[grepl("Populus_deltoides", phy.sps.uniqu)]
matchednamessilvics$sppMatch[which(matchednamessilvics$silvicsname == "Populus_deltoides_subsp._monilifera")] <-  "Populus_deltoides"

#phy.sps.uniqu[grepl("Acer_saccharum", phy.sps.uniqu)]
matchednamessilvics$sppMatch[which(matchednamessilvics$silvicsname == "Acer_glabrum_var._glabrum")] <-  "Acer_saccharum_subsp._saccharum"

phy.sps.uniqu[grepl("Schefflera_morototoni", phy.sps.uniqu)]
matchednamessilvics$sppMatch[which(matchednamessilvics$silvicsname == "Didymopanax_morototoni")] <-  "Schefflera_morototoni"

#phy.sps.uniqu[grepl("Eucalyptus_globulus", phy.sps.uniqu)]
matchednamessilvics$sppMatch[which(matchednamessilvics$silvicsname == "Eucalyptus_globulus")] <-  "Eucalyptus_globulus_subsp._globulus"

#phy.sps.uniqu[grepl("Magnolia_virginiana", phy.sps.uniqu)]
matchednamessilvics$sppMatch[which(matchednamessilvics$silvicsname == "Magnolia_virginiana")] <-  "Magnolia_virginiana_var._australis"

#phy.sps.uniqu[grepl("", phy.sps.uniqu)]
matchednamessilvics$sppMatch[which(matchednamessilvics$silvicsname == "Taxodium_distichum_var._distichum")] <-  "Taxodium_distichum"

nomatchAfterKewcheck <- subset(matchednamessilvics , is.na(sppMatch))
nrow(nomatchAfterKewcheck)

matchnona <- matchednamessilvics[!is.na(matchednamessilvics$sppMatch),]
d$sppMatch <- NA
# first add species that match the tree already
easymatch <- silvics.sps[which(silvics.sps %in% phy.sps.uniqu)] 
d$sppMatch <- ifelse(d$latbi %in% easymatch, d$latbi, NA)

missing_idx <- is.na(d$sppMatch)
matched_values <- match(d$latbi[missing_idx], matchnona$silvicsname)

# match spp names back in usda df
d$sppMatch[missing_idx] <- matchnona$sppMatch[matched_values]
# check
d[!duplicated(d$latbi), c("latbi", "sppMatch")]

silvics_sub <- d[c("latbi","genusName","sppMatch")]

set.seed(123)
vec <- "Tilia"
t <- subset(silvics_sub, genusName %in% vec)
t <- t$latbi
spp_smalltree <- unique(t)
studyNosmall <- subset(silvics_sub, latbi %in% spp_smalltree)


smallTree <- keep.tip(phy.plants, which(phy.plants$tip.label %in% unique(t)))

smallnamesphy <- smallTree$tip.label
smallTree$root.edge <- 0

is.rooted(smallTree)
smallTree$node.label<-NULL

studyNosmall$latbi[which(studyNosmall$latbi %in% smallTree$tip.label)]

smallTreeUltra <- force.ultrametric(
  smallTree,
  method = "extend"  # Extends terminal branches
)

unique(nomatch$silvicsname)
spptoplice <- c()

for (genus in vec) {
  matches <- nomatchAfterKewcheck$silvicsname[grepl(genus, nomatchAfterKewcheck$silvicsname)]
  spptoplice <- c(spptoplice, matches)
}

smallTreeUltraSpliced <- smallTreeUltra

# Loop and update the tree with one species a time
for (i in spptoplice) {
  smallTreeUltraSpliced <- add.species.to.genus(tree = smallTreeUltraSpliced, species = i, where = "root")
}

smallnamesphy <- smallTreeUltraSpliced$tip.label
studyNosmallSpliced <- subset(d, latbi %in% smallnamesphy)

smallTreeUltraSpliced$root.edge <- 0

is.rooted(smallTreeUltraSpliced)
smallTreeUltraSpliced$node.label<-NULL

smallTreeUltraSpliced
studyNosmallSpliced

# remove duplicated rows
studyNosmallSpliced <- studyNosmallSpliced[!duplicated(studyNosmallSpliced$latbi), ]

rownames(studyNosmallSpliced) <- studyNosmallSpliced$latbi

# species in both tree and data
common_species <- intersect(smallTreeUltraSpliced$tip.label, studyNosmallSpliced$latbi)

# drop unmatched species from the tree and data
pruned_tree <- drop.tip(smallTreeUltraSpliced, setdiff(smallTreeUltraSpliced$tip.label, common_species))

# list of species got spliced
target <- nomatch$silvicsname

# get the index of the tip
tip_index <- match(target, pruned_tree$tip.label)
# make the full usda tree
silvicslist <- unique(d$sppMatch)

silvicsTree <- keep.tip(phy.plants, phy.plants$tip.label[phy.plants$tip.label %in% silvicslist])
matchednamessilvics1 <- matchednamessilvics[!is.na(matchednamessilvics$sppMatch), ]
name_map <- setNames(matchednamessilvics1$silvicsname, matchednamessilvics1$sppMatch)

# replace the tip name with the name in usda
silvicsTree$tip.label <- ifelse(silvicsTree$tip.label %in% names(name_map),
                                name_map[silvicsTree$tip.label],
                                silvicsTree$tip.label)

silvicsUltra <- force.ultrametric(
  silvicsTree,
  method = "extend"  # Extends terminal branches
)

# splicing
silvicsSpliced <- silvicsUltra
for (i in spptoplice) {
  silvicsSpliced <- add.species.to.genus(tree = silvicsSpliced, species = i, where = "root")
}

silvicsSpliced$tip.label[duplicated(silvicsSpliced$tip.label)]
unique(silvicsSpliced$tip.label)
unique(d$latbi)

setdiff(d$latbi,silvicsSpliced$tip.label)
# write out the tree
write.tree(silvicsSpliced,"output/silvicsPhylogenyFull.tre")
}

silvicsTree <- read.tree("output/silvicsPhylogenyFull.tre")
masting <- d[!duplicated(d$latbi), c("latbi","mastEvent")]
masting$mastEvent[is.na(masting$mastEvent)] <- "No information"

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

preliminaryPhylo <- FALSE
if(preliminaryPhylo){
# Assign colors to each dormancy class
mastingYN <- unique(masting$mastEvent)
mast_colors <- c("Y" = "#CD5555", "N" = "#698B22", "No information" = "Grey")

masting$color <- mast_colors[masting$mastEvent]

pdf("output/figures/mastTraitTree.pdf", width = 20, height = 20)
plot(
  silvicsTree, type = 'fan', label.offset = 0.2,
  cex = 0.6, no.margin = TRUE, tip.color = masting$color
)


# Add monoecious or dioecious
monodio <- subset(d, !is.na(d$typeMonoOrDio))
monodiosp <- monodio[!duplicated(monodio$latbi), c("latbi","typeMonoOrDio")]

# Match species to tips
tipindex_monodio <- match(monodiosp$latbi, allspp)

monocolor <- c("Monoecious" = "#659794", "Dioecious" = "#A6739B", "Polygamous" = "#D4C95F")

# Add symbols at tips
for (i in seq_along(monodiosp$latbi)) {
  sp <- monodiosp$latbi[i]
  tip <- tipindex_monodio[i]
  type <- monodiosp$typeMonoOrDio[i]
  # Get color and shape
  col <- monocolor[type]
  pch <- 21
  tiplabels(pch = pch, tip = tip,  offset=75, col = "black", bg = col, cex = 1.2)
}

# Add seed dormancy
dormancy <- subset(d, !is.na(d$seedDormancy))
dormancysp <- dormancy[!duplicated(dormancy$latbi), c("latbi","seedDormancy")]

# Match species to tips
tipindex_dormancy <- match(dormancysp$latbi, allspp)

dormcolor <- c("Y" = "#5BA2CC", "N" = "#F2BA68")

# Add symbols at tips
for (i in seq_along(dormancysp$latbi)) {
  sp <- dormancysp$latbi[i]
  tip <- tipindex_dormancy[i]
  type <- dormancysp$seedDormancy[i]
  # Get color and shape
  col <- dormcolor[type]
  pch <- 21
  tiplabels(pch = pch, tip = tip,  offset=85, col = "black", bg = col, cex = 1.2)
}

# Add pollination mode
pollination <- subset(d, !is.na(d$pollination))
pollinationsp <- pollination[!duplicated(pollination$latbi), c("latbi","pollination")]

pollination_list <- strsplit(as.character(pollinationsp$pollination), "and\\s*")
species <- as.character(pollinationsp$latbi)
# Match species to tips
tipindex_pollination <- match(species, allspp)

pollinationcolor <- c("wind" = "#AC7299", "animals" = "#CFDAA8")

base_offset <- 95         # where symbols begin
symbol_spacing <- 5       # space between multiple symbols
symbol_size <- 1.2        # size of symbols
# Add symbols at tips
for (i in seq_along(species)) {
  tip <- tipindex_pollination[i]
  modes <- pollination_list[[i]]
  n_modes <- length(modes)
  
  # Calculate centered offsets
  offsets <- base_offset + ((seq_len(n_modes) - (n_modes + 1)/2) * symbol_spacing)
  
  for (j in seq_along(modes)) {
    mode <- modes[j]
    col <- pollinationcolor[mode]
    
    # Plot symbol
    tiplabels(pch = pch, tip = tip, offset = offsets[j], col = "black", bg = col, cex = symbol_size)
  }
}

# Add seed dispersal mode
dispersal <- subset(d, !is.na(d$seedDispersal))
dispersalsp <- dispersal[!duplicated(dispersal$latbi), c("latbi","seedDispersal")]

dispersal_list <- strsplit(as.character(dispersalsp$seedDispersal), ",\\s*")
species <- as.character(dispersalsp$latbi)
# Match species to tips
tipindex_dispersal <- match(species, allspp)

dispersalcolor <- c("wind" = "#C43142","gravity" = "#C43142", "water" ="#C43142", "fire" = "#C43142", "mammals" = "#387E46", "birds" = "#387E46", "rodents" = "#387E46")

base_offset <- 110         # where symbols begin
symbol_spacing <- 5       # space between multiple symbols
symbol_size <- 1.2        # size of symbols
# Add symbols at tips
for (i in seq_along(species)) {
  tip <- tipindex_dispersal[i]
  modes <- dispersal_list[[i]]
  n_modes <- length(modes)
  
  # Calculate centered offsets
  offsets <- base_offset + ((seq_len(n_modes) - (n_modes + 1)/2) * symbol_spacing)
  
  for (j in seq_along(modes)) {
    mode <- modes[j]
    col <- dispersalcolor[mode]
    
    # Plot symbol
    tiplabels(pch = pch, tip = tip, offset = offsets[j], col = "black", bg = col, cex = symbol_size)
  }
}

# Add seed Weight
seedWeight <- subset(d, !is.na(d$seedWeights))
seedWeightsp <- monodio[!duplicated(seedWeight$latbi), c("latbi","seedWeights")]
values <- na.omit(log(seedWeightsp$seedWeight))

# Normalize the values
scaled_cex <- scales::rescale(values, to = c(0.5, 4)) 

# Match species to tips
tipindex_seedweight <- match(seedWeight$latbi, allspp)

# Plot variable-size symbols
for (i in seq_along(seedWeight$latbi)) {
  tip <- tipindex_seedweight[i]
  size <- scaled_cex[i]
  
  tiplabels(
    pch = 21,                     # Circle shape
    tip = tip, 
    offset = 130,
    col = "black", 
    bg = "#E87687",              # Or another color if you like
    cex = size
  )
}

# Add legend
legend(
  x=-90, y=100,
  legend = names(mast_colors),
  col = mast_colors,
  pch = 15,           
  title = "Masting",
  pt.cex = 1.2,          # size of legend symbols
  cex = 1,             # text size
  bty = "n"            # no box around legend
)
# Reproductive system legend
legend(x=70, y=40,
       legend = names(monocolor),
       pt.bg = unname(monocolor), 
       pch = 21,        
       pt.cex = 1.2,               
       cex = 1,                  
       title = "Reproductive system",
       bty = "n")
# Seed dormancy legend
legend(x=70, y=100,
       legend = names(dormcolor),
       pt.bg = unname(dormcolor), 
       pch = 21,        
       pt.cex = 1.2,               
       cex = 1,                  
       title = "Seed Dormancy",
       bty = "n")
# pollination mode legend
legend(x=70, y=-10, 
       legend = names(pollinationcolor),
       pt.bg = pollinationcolor, 
       pch = 21,
       pt.cex = 1.2,
       cex = 1,
       title = "Pollination Mode",
       bty = "n")  

# dispersal mode legend
legend(x=70, y=-50, 
       legend = names(dispersalcolor),
       pt.bg = dispersalcolor, 
       pch = 21,
       pt.cex = 1.2,
       cex = 1,
       title = "Dispersal Mode",
       bty = "n") 
# Seed weight legend
legend_sizes <- c(min(values), median(values), max(values))
legend_cex <- scales::rescale(legend_sizes, to = c(0.5, 2))

legend(
  x = -90, y = 40,
  legend = round(legend_sizes, 2),
  pt.cex = legend_cex,
  pch = 21,
  pt.bg = "#E87687",
  col = "black",
  title = "Fruit size",
  bty = "n"
)

dev.off()
}
quercus <- FALSE
if(quercus){# Quercus tree
phy.plants<-read.tree("data/quercusFullTree.tre")
phy.plants$tip.label <- sapply(phy.plants$tip.label, function(x) {
  parts <- strsplit(x, "[_|]")[[1]] 
  paste(parts[1], parts[2], sep = "_")
})
phy.plants$tip.label
quercus <- grep("^Quercus", phy.plants$tip.label, value = TRUE)
quercusTree <- keep.tip(phy.plants, quercus)

silvicsQuercus <- d[grep("^Quercus", d$latbi), ]
quercusY <- silvicsQuercus$latbi[silvicsQuercus$mastEvent == "Y"]
quercusN  <- silvicsQuercus$latbi[silvicsQuercus$mastEvent == "N"]
tipColors <- rep("black", length(quercusTree$tip.label))
tipColors[quercusTree$tip.label %in% quercusY] <- "#ED562C"
tipColors[quercusTree$tip.label %in% quercusN] <- "#A9D5B1"

pdf("output/figures/quercusTree.pdf", width = 100, height = 100)

plot(
  quercusTree, ,type = "fan",label.offset = 0.001,
  cex = 1.2, no.margin = TRUE, tip.color = tipColors
)

legend("bottomright",
       legend = c("Masting", "Non-masting"),
       col = c("#ED562C", "#A9D5B1"),
       pch = 19,
       bty = "n",
       cex = 5)

dev.off()
}
# Calculating phylogenetic signal -- angiosperm
# Continuous
weight <- angio$logSeedWeight
names(weight) <- angio$latbi
weightLambda <- phylosig(phyangio, weight, method = "lambda", test = TRUE)
print(weightLambda)
# lambda is almost 1
fruit <- angio$logFruit
names(fruit) <- angio$latbi
fruitLambda <- phylosig(phyangio, fruit, method = "lambda", test = TRUE)
print(fruitLambda)
# lambda is almost 1
seed <- angio$logSeedSize
names(seed) <- angio$latbi
seedLambda <- phylosig(phyangio, seed, method = "lambda", test = TRUE)
print(seedLambda)

oil <- angio$oilContent
names(oil) <- angio$latbi
oilLambda <- phylosig(phyangio, oil, method = "lambda", test = TRUE)
print(oilLambda)

leaf <- angio$leafLongevity
names(leaf) <- angio$latbi
leafLambda <- phylosig(phyangio, leaf, method = "lambda", test = TRUE)
print(leafLambda)
# lambda is almost 1

# categorical
phyangio <- multi2di(phyangio)
disp <- angio$seedDispersal
names(disp) <- angio$latbi
dispLambda <- fitDiscrete(phyangio, disp, model="ER",transform="lambda")
print(dispLambda$opt$lambda)

# lambda is 1
dorm <- angio$seedDormancy
names(dorm) <- angio$latbi
dormLambda <- fitDiscrete(phyangio, dorm, model="ER",transform="lambda")


# lambda is close to 0
poll <- angio$pollination
names(poll) <- angio$latbi
pollLambda <- fitDiscrete(phyangio, poll, model="ER",transform="lambda")
# lambda is close to 1
rep <- angio$typeMonoOrDio
names(rep) <- angio$latbi
repLambda <- fitDiscrete(phyangio, rep, model="ER",transform="lambda")
# lambda is close to 0
drought <- angio$droughtTolerance
names(drought) <- angio$latbi
droughtLambda <- fitDiscrete(phyangio, drought, model="ER",transform="lambda")
# lambda is close to 0
mast <- angio$mastEvent
names(mast) <- angio$latbi
mastLambda <- fitDiscrete(phyangio, mast, model="ER",transform="lambda")
print(mastLambda$opt$lambda)
# D statistics
set.seed(135)
phyangio   <- drop.tip(silvicsTree, setdiff(silvicsTree$tip.label, angio$latbi))
phyangio$node.label <- NULL
dMast <- angio[, c("mastEvent", "latbi")]
comp <- comparative.data(phyangio, dMast, names.col="latbi", vcv=TRUE)
phylo.d(comp, binvar=mastEvent)
#Estimated D :  0.5779433
#Probability of E(D) resulting from no (random) phylogenetic structure :  0
#Probability of E(D) resulting from Brownian phylogenetic structure    :  0.01

dDorm <- angio[, c("seedDormancy", "latbi")]
comp <- comparative.data(phyangio, dDorm, names.col="latbi", vcv=TRUE)

phylo.d(comp, binvar=seedDormancy)
#Estimated D :  0.2540649
#Probability of E(D) resulting from no (random) phylogenetic structure :  0
#Probability of E(D) resulting from Brownian phylogenetic structure    :  0.161

lambda_angio <- data.frame(
  Group = 
    rep("Angiosperm", 8)
  ,
  Trait = c("Seed weight (Log)", "Fruit size (Log)",
    "Oil content", "Leaf longevity",
    "Dispersal mode", 
    "Pollination mode", "Reproductive type", "Drought tolerance"),
  Lambda = c(
    round(weightLambda$lambda,3),
    round(fruitLambda$lambda,3),
    round(oilLambda$lambda,3),
    round(leafLambda$lambda,3),
    round(dispLambda$opt$lambda,3),
    round(pollLambda$opt$lambda,3),
    round(repLambda$opt$lambda,3),
    round(droughtLambda$opt$lambda,3)
  ),
  stringsAsFactors = FALSE
)

d_angio <- data.frame(
  Group = 
    rep("Angiosperm", 2)
  ,
  Trait = c("Seed dormancy", "Masting"),
  "Estimated D" = c("0.25","0.58"
  ),"P(random)" = c("0","0"),"P(Brownian)" = c("0.16","0.01"),
  
  stringsAsFactors = FALSE,check.names = FALSE
)

# Predation satiation for angiosperm, using seed weight for the branch color

commonSp <- intersect(phyangio$tip.label, names(weight))
mast <- angio$mastEvent
names(mast) <- angio$latbi
mast <- mast[!is.na(mast)]
commonSp <- intersect(commonSp,names(mast))
weight <- weight[commonSp]
weight <- weight[!is.na(weight)]
weightTree <- drop.tip(phyangio,
                        setdiff(phyangio$tip.label,
                                names(weight)))
name.check(weightTree, weight)



pdf("output/figures/weightTree1.pdf", width = 10, height = 20)
weightMap <- contMap(weightTree,
                     weight,fsize=c(1,1),
                     plot = TRUE)
lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
## for fun, let's change our contMap gradient
object <- setMap(weightMap, colors = colorRampPalette(c("navy", "white", "darkred"))(100))
plot(object,ftype="off",fzie=c(0.5,0.7),xlim = c(0, max(lastPP$xx) * 6),
     leg.txt="Seed weight (Log)",legend=100,sig=1,lwd=2)
lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
## make sure our discrete character is in the order of our tree
weightDf <- d[!is.na(d$seedWeights), ]
weightMast <- weightDf$mastEvent
names(weightMast) <- weightDf$latbi
weightMast  <- weightMast[weightTree$tip.label]
cols <- setNames(
  c("#A9D5B1", "#ED562C"),
  levels(factor(weightMast))
)
weightDisp <- weightDf$seedDispersal
names(weightDisp) <- weightDf$latbi
weightDisp  <- weightDisp[weightTree$tip.label]
colsDisp <- setNames(
  c("#009392", "#72AAA1", "#B1C7B3",na.value = "white"),
  levels(factor(weightDisp))
)
weightOil <- weightDf$oilContent
names(weightOil) <- weightDf$latbi
weightOil  <- weightOil[weightTree$tip.label]
colsOil <- "#FCD116"

weightDorm <- weightDf$seedDormancy
names(weightDorm) <- weightDf$latbi
weightDorm  <- weightDorm[weightTree$tip.label]
colsDorm <- setNames(
  c("#CFDAA8","#7F883B", na.value = "white"),
  levels(factor(weightDorm))
)

weightRep <- weightDf$typeMonoOrDio
names(weightRep) <- weightDf$latbi
weightRep  <- weightRep[weightTree$tip.label]
colsRep <- setNames(
  c("#D0587E", "#D98994", "#E5B9AD",na.value = "white"),
  levels(factor(weightRep))
)

weightPoll <- weightDf$pollination
names(weightPoll) <- weightDf$latbi
weightPoll  <- weightPoll[weightTree$tip.label]
colsPoll <- setNames(
  c("Animal" = "#08519C","Both" = "#4292C6", "Wind" = "#92d4f5",na.value = "white"),
  levels(factor(weightPoll))
)

weightDrought <- weightDf$droughtTolerance
names(weightDrought) <- weightDf$latbi
weightDrought <- weightDrought[weightTree$tip.label]
colsDrought <- setNames(
  c("High" = "#da691f","Moderate" = "#e79259", "Low" = "#f0bb98",na.value = "white"),
  levels(factor(weightDrought))
)

weightLeaf <- weightDf$leafLongevity
names(weightLeaf) <- weightDf$latbi
weightLeaf  <- weightLeaf[weightTree$tip.label]
colsLeaf <- "#5E3C99"

weightMast <- weightMast[object$tree$tip.label]
weightDisp <- weightDisp[object$tree$tip.label]
weightOil <- weightOil[object$tree$tip.label]
weightDorm <- weightDorm[object$tree$tip.label]
weightRep <- weightRep[object$tree$tip.label]
weightPoll <- weightPoll[object$tree$tip.label]
weightDrought <- weightDrought[object$tree$tip.label]
weightLeaf <- weightLeaf[object$tree$tip.label]

Ntip <- length(object$tree$tip.label)
tree_width <- diff(range(lastPP$xx))
for(i in 1:length(weightMast)){
  text(lastPP$xx[i],
       lastPP$yy[i],
       gsub("_", " ", weightMap$tree$tip.label[i]),
       pos=4,
       col=cols[weightMast[i]],
       cex=0.7,
       font=3)
}

for(i in 1:length(weightDisp)){
points(lastPP$xx[1:Ntip] + tree_width * 1.3,
       lastPP$yy[1:Ntip],
       pch = 21,
       bg = colsDisp[weightDisp],
       col = "white",
       cex = 1.2)
}

sizeMin <- 0.6
sizeMax <- 1.2

oilSize <- (weightOil - min(weightOil, na.rm=TRUE)) /
  (max(weightOil, na.rm=TRUE) - min(weightOil, na.rm=TRUE))

oilSize <- sizeMin + oilSize * (sizeMax - sizeMin)
  for(i in 1:length(weightOil)){
  points(lastPP$xx[1:Ntip] + tree_width *1.6,
         lastPP$yy[1:Ntip],
         pch=21,
         bg="#FCD116",
         col="white",
         cex=oilSize)
}

for(i in 1:length(weightDorm)){
  points(lastPP$xx[1:Ntip] + tree_width *1.9,
         lastPP$yy[1:Ntip],
         pch=21,
         bg=colsDorm[weightDorm],
         col="white",
         cex=1.2)
}

for(i in 1:length(weightRep)){
  points(lastPP$xx[1:Ntip] + tree_width *2.3,
         lastPP$yy[1:Ntip],
         pch=21,
         bg=colsRep[weightRep],
         col="white",
         cex=1.2)
}

for(i in 1:length(weightPoll)){
  points(lastPP$xx[1:Ntip] + tree_width *2.6,
         lastPP$yy[1:Ntip],
         pch=21,
         bg=colsPoll[weightPoll],
         col="white",
         cex=1.2)
}

for(i in 1:length(weightDrought)){
  points(lastPP$xx[1:Ntip] + tree_width *3,
         lastPP$yy[1:Ntip],
         pch=21,
         bg=colsDrought[weightDrought],
         col="white",
         cex=1.2)
}

sizeMin <- 0.6
sizeMax <- 1.2

leafSize <- (weightLeaf - min(weightLeaf, na.rm=TRUE)) /
  (max(weightLeaf, na.rm=TRUE) - min(weightLeaf, na.rm=TRUE))

leafSize <- sizeMin + leafSize * (sizeMax - sizeMin)
for(i in 1:length(weightLeaf)){
  points(lastPP$xx[1:Ntip] + tree_width *3.3,
         lastPP$yy[1:Ntip],
         pch=21,
         bg="#5E3C99",
         col="white",
         cex=leafSize)
}
usr <- par("usr")

# Title
text(x = usr[2] * 0.21,
     y = usr[4] * 0.98, "Masting", pos = 2,cex = 0.7)

# Legend right next to it
legend(x = usr[2] * 0.23,
       y = usr[4] * 0.98,
       legend = c("Y", "N"),
       fill = cols,
       border = NA,
       bty = "n",
       horiz = FALSE,
       xjust = 0,
       yjust = 0.5,cex = 0.7)

text(x = usr[2] * 0.37, y = usr[4] * 0.99, labels = paste(c("Dispersal", "mode"), collapse = "\n"), adj = c(0.5,0.5), cex = 0.3)

legend(x = usr[2] * 0.35,
       y = usr[4] * 0.975,
       legend = c("Abiotic", "Biotic", "Both"),
       fill = c("#009392", "#72AAA1", "#B1C7B3"),
       border = NA,
       bty = "n",
       horiz = FALSE,
       xjust = 0,
       yjust = 0.5,cex = 0.3)

legend_vals <- c(min(weightOil, na.rm=TRUE),
                 mean(weightOil, na.rm=TRUE),
                 max(weightOil, na.rm=TRUE))
legend_sizes <- (legend_vals - min(weightOil, na.rm=TRUE)) /
  (max(weightOil, na.rm=TRUE) - min(weightOil, na.rm=TRUE))

legend_sizes <- sizeMin + legend_sizes * (sizeMax - sizeMin)

text(x = usr[2] * 0.42, y = usr[4] * 0.99, labels = paste(c("Oil", "content"), collapse = "\n"), adj = c(0.5,0.5), cex = 0.3)
legend(x = usr[2] * 0.4,
       y = usr[4] * 0.983,
       legend = round(legend_vals, 2),
       pt.cex = legend_sizes,
       pch = 21,
       pt.bg = "#FCD116",
       col = "white",
       bty = "n",
       horiz = FALSE,
       cex = 0.3)

text(x = usr[2] * 0.47, y = usr[4] * 0.99, labels = paste(c("Seed", "dormancy"), collapse = "\n"), adj = c(0.5,0.5), cex = 0.3)

legend(x = usr[2] * 0.44,
       y = usr[4] * 0.975,
       legend = c("Dormant", "Non-dormant"),
       fill = c("#7F883B", "#CFDAA8"),
       border = NA,
       bty = "n",
       horiz = FALSE,
       xjust = 0,
       yjust = 0.5,cex = 0.3)



text(x = usr[2] * 0.53, y = usr[4] * 0.99, labels = paste(c("Reproductive", "mode"), collapse = "\n"), adj = c(0.5,0.5), cex = 0.3)


legend(x = usr[2] * 0.5,
       y = usr[4] * 0.975,
       legend = c("Dioecious", "Monoecious", "Polygamous"),
       fill = c("#D0587E", "#D98994","#E5B9AD"),
       border = NA,
       bty = "n",
       horiz = FALSE,
       xjust = 0,
       yjust = 0.5,cex = 0.3)

text(x = usr[2] * 0.58, y = usr[4] * 0.99, labels = paste(c("Pollination", "mode"), collapse = "\n"), adj = c(0.5,0.5), cex = 0.3)

legend(x = usr[2] * 0.56,
       y = usr[4] * 0.975,
       legend = c("Animal", "Both", "Wind"),
       fill = c("#08519C", "#4292C6","#92d4f5"),
       border = NA,
       bty = "n",
       horiz = FALSE,
       xjust = 0,
       yjust = 0.5,cex = 0.3)

text(x = usr[2] * 0.64, y = usr[4] * 0.99, labels = paste(c("Drought", "tolerance"), collapse = "\n"), adj = c(0.5,0.5), cex = 0.3)

legend(x = usr[2] * 0.62,
       y = usr[4] * 0.975,
       legend = c("High", "Moderate", "Low"),
       fill = c("#da691f", "#e79259","#f0bb98"),
       border = NA,
       bty = "n",
       horiz = FALSE,
       xjust = 0,
       yjust = 0.5,cex = 0.3)

legend_vals <- c(min(weightLeaf, na.rm=TRUE),
                 mean(weightLeaf, na.rm=TRUE),
                 max(weightLeaf, na.rm=TRUE))
legend_sizes <- (legend_vals - min(weightLeaf, na.rm=TRUE)) /
  (max(weightLeaf, na.rm=TRUE) - min(weightLeaf, na.rm=TRUE))

legend_sizes <- sizeMin + legend_sizes * (sizeMax - sizeMin)

text(x = usr[2] * 0.69, y = usr[4] * 0.99, labels = paste(c("Leaf", "longevity"), collapse = "\n"), adj = c(0.5,0.5), cex = 0.3)
legend(x = usr[2] * 0.68,
       y = usr[4] * 0.983,
       legend = round(legend_vals, 2),
       pt.cex = legend_sizes,
       pch = 21,
       pt.bg = "#5E3C99",
       col = "white",
       bty = "n",
       horiz = FALSE,
       cex = 0.3)


text(x = usr[2] * 0.85,
     y = usr[4] * 0.01,
     labels = paste(c("Lambda:",
                      "Seed weight (Log): 0.95",
                      "Fruit size (Log): 0.93",
                      "Dispersal mode: 1",
                      "Oil content: 0.25",
                      "Leaf longevity: 1.02",
                      "Pollination mode: 1",
                      "Reproductive type: 1",
                      "Drought tolerance: 0",
                      "D-statistic", 
                      "Seed dormancy: 0.25; P(random) = 0; P(Brownian) = 0.16",
                      "Masting: 0.57; P(random) = 0; P(Brownian) = 0.01"),
                    collapse = "\n"),
     pos = 3, cex =0.7)

#text(x = 285, y = 92,
#     labels = paste(c("Predation",
#                      "satiation"),collapse = "\n"),
#     pos = 3, cex =0.7)

#text(x = 345, y = 92,
#     labels = paste(c("Pollination",
#                      "coupling"),collapse = "\n"),
 #    pos = 3, cex =0.7)

#text(x = 405, y = 92,
#     labels = paste(c("Resource",
#                      "matching"),collapse = "\n"),
#     pos = 3, cex =0.7)

dev.off()

# Calculating phylogenetic signal -- gymnosperm
weight <- conifer$logSeedWeight
names(weight) <- conifer$latbi
weightLambda <- phylosig(phyconifer, weight, method = "lambda", test = TRUE)
print(weightLambda)
# lambda is almost 1
fruit <- conifer$logFruit
names(fruit) <- conifer$latbi
fruitLambda <- phylosig(phyconifer, fruit, method = "lambda", test = TRUE)
print(fruitLambda)
# lambda is almost 1
seed <- conifer$logSeedSize
names(seed) <- conifer$latbi
seedLambda <- phylosig(phyconifer, seed, method = "lambda", test = TRUE)
print(seedLambda)

oil <- conifer$oilContent
names(oil) <- conifer$latbi
oilLambda <- phylosig(phyconifer, oil, method = "lambda", test = TRUE)
print(oilLambda)
# lambda is close to 0
leaf <- conifer$leafLongevity
names(leaf) <- conifer$latbi
leafLambda <- phylosig(phyconifer, leaf, method = "lambda", test = TRUE)
print(leafLambda)
# lambda is 0.896

phyconifer <- multi2di(phyconifer)
disp <- conifer$seedDispersal
names(disp) <- conifer$latbi
dispLambda <- fitDiscrete(phyconifer, disp, model="ER",transform="lambda")
print(dispLambda$opt$lambda)
# lambda is 0
dorm <- conifer$seedDormancy
names(dorm) <- conifer$latbi
dormLambda <- fitDiscrete(phyconifer, dorm, model="ER",transform="lambda")
print(dormLambda$opt$lambda)
# lambda is close to 0

rep <- conifer$typeMonoOrDio
names(rep) <- conifer$latbi
repLambda <- fitDiscrete(phyconifer, rep, model="ER",transform="lambda")
print(repLambda$opt$lambda)
# lambda is close to 1
drought <- conifer$droughtTolerance
names(drought) <- conifer$latbi
droughtLambda <- fitDiscrete(phyconifer, drought, model="ER",transform="lambda")
print(droughtLambda$opt$lambda)
# lambda is close to 1
mast <- conifer$mastEvent
names(mast) <- conifer$latbi
mastLambda <- fitDiscrete(phyconifer, mast, model="ER",transform="lambda")
print(droughtLambda$opt$lambda)

# D statistics
phyconifer <- drop.tip(silvicsTree, setdiff(silvicsTree$tip.label, conifer$latbi))
phyconifer$node.label<-NULL
dDorm <- conifer[, c("seedDormancy", "latbi")]
comp <- comparative.data(phyconifer, dDorm, names.col="latbi", vcv=TRUE)
phylo.d(comp, binvar=seedDormancy)
# Estimated D :  0.9764245
# Probability of E(D) resulting from no (random) phylogenetic structure :  0.412
# Probability of E(D) resulting from Brownian phylogenetic structure    :  0
dRep <- conifer[, c("typeMonoOrDio", "latbi")]
comp <- comparative.data(phyconifer, dRep, names.col="latbi", vcv=TRUE)
phylo.d(comp, binvar=typeMonoOrDio)
# Estimated D :  -0.5188342
# Probability of E(D) resulting from no (random) phylogenetic structure :  0
# Probability of E(D) resulting from Brownian phylogenetic structure    :  0.885
dMast <- conifer[, c("mastEvent", "latbi")]
comp <- comparative.data(phyconifer, dMast, names.col="latbi", vcv=TRUE)
phylo.d(comp, binvar=mastEvent)
#Estimated D :  0.5059193
#Probability of E(D) resulting from no (random) phylogenetic structure :  0.02
#Probability of E(D) resulting from Brownian phylogenetic structure    :  0.097
lambda_conifer <- data.frame(
  Group = 
    rep("Gymnosperm", 6)
  ,
  Trait = c("Seed weight (Log)", "Fruit size (Log)",
            "Oil content", "Leaf longevity",
            "Dispersal mode", 
            "Drought tolerance"),
  Lambda = c(
    round(weightLambda$lambda,3),
    round(fruitLambda$lambda,3),
    round(oilLambda$lambda,3),
    round(leafLambda$lambda,3),
    round(dispLambda$opt$lambda,3),
    round(droughtLambda$opt$lambda,3)
  ),
  stringsAsFactors = FALSE
)
d_conifer <- data.frame(
  Group = 
    rep("Gymnosperm", 3)
  ,
  Trait = c("Masting", "Seed dormancy", "Reproductive type"),
  "Estimated D" = c("0.51","0.98", "-0.52"
  ),"P(random)" = c("0.02","0.41", "0"),"P(Brownian)" = c("0.10","0", "0.89"),
  
  stringsAsFactors = FALSE,check.names = FALSE
)

lambda <- rbind(lambda_angio,lambda_conifer)
dStat <- rbind(d_angio,d_conifer)
labmdaTable <- xtable(lambda, 
                      caption = "Phylogenetic Signal (Pagel's ${lambda}$) for continous traits and multi-level categorical traits", 
                      label = "tab:lambda")
print(labmdaTable, type = "latex", include.rownames = FALSE)
dTable <- xtable(dStat, 
                 caption = "Phylogenetic Signal (D-statistic with probability from random phylogenetic structure and probability from Brownian phylogenetic structure for binary traits", 
                 label = "tab:d")
print(dTable, type = "latex", include.rownames = FALSE)