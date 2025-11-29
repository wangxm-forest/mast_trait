#### Phylogeny####
####Start by Mao####
####Sep 17####

library(ape)
library(phytools)
library(viridisLite)
# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

setwd("C:/PhD/Project/PhD_thesis/mast_trait/")
source("analyses/dataCleaning.R")

maketree <- FALSE
if(maketree){
phy.plants<-read.tree("C:/PhD/Project/egret/analyses/input/ALLMB.tre")

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

silvics.phenosp.sps.inphylo <- silvics.sps[which(!silvics.sps%in%phy.sps.uniqu)] #20 out of 190

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

phy.sps.uniqu[grepl("Notholithocarpus_densiflorus", phy.sps.uniqu)]
matchednamessilvics$sppMatch[which(matchednamessilvics$silvicsname == "Notholithocarpus_densiflorus")] <-  "Notholithocarpus_densiflorus_var._densiflorus"

#phy.sps.uniqu[grepl("Calophyllum_calaba", phy.sps.uniqu)]
matchednamessilvics$sppMatch[which(matchednamessilvics$silvicsname == "Calophyllum_calaba")] <-  "Calophyllum_calaba_var._bracteatum"

#phy.sps.uniqu[grepl("Populus_deltoides", phy.sps.uniqu)]
matchednamessilvics$sppMatch[which(matchednamessilvics$silvicsname == "Populus_deltoides_subsp._monilifera")] <-  "Populus_deltoides"

phy.sps.uniqu[grepl("Schefflera_morototoni", phy.sps.uniqu)]
matchednamessilvics$sppMatch[which(matchednamessilvics$silvicsname == "Didymopanax_morototoni")] <-  "Schefflera_morototoni"

#phy.sps.uniqu[grepl("Acer", phy.sps.uniqu)]
matchednamessilvics$sppMatch[which(matchednamessilvics$silvicsname == "Acer_nigrum")] <-  "Acer_saccharum_subsp._nigrum"

#phy.sps.uniqu[grepl("Eucalyptus_globulus", phy.sps.uniqu)]
matchednamessilvics$sppMatch[which(matchednamessilvics$silvicsname == "Eucalyptus_globulus")] <-  "Eucalyptus_globulus_subsp._globulus"

#phy.sps.uniqu[grepl("Magnolia_virginiana", phy.sps.uniqu)]
matchednamessilvics$sppMatch[which(matchednamessilvics$silvicsname == "Magnolia_virginiana")] <-  "Magnolia_virginiana_var._australis"

#phy.sps.uniqu[grepl("Quercus_falcata_var._pagodifolia", phy.sps.uniqu)]
matchednamessilvics$sppMatch[which(matchednamessilvics$silvicsname == "Quercus_falcata_var._pagodifolia")] <-  "Quercus_pagoda"

#phy.sps.uniqu[grepl("", phy.sps.uniqu)]
matchednamessilvics$sppMatch[which(matchednamessilvics$silvicsname == "Taxodium_distichum_var._distichum")] <-  "Taxodium_distichum"

#phy.sps.uniqu[grepl("Quercus_falcata", phy.sps.uniqu)]
matchednamessilvics$sppMatch[which(matchednamessilvics$silvicsname == "Quercus_falcata_var._falcata")] <-  "Quercus_falcata"

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

allspp <- silvicsTree$tip.label

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

