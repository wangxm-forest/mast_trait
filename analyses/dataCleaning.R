### Data cleaning ###
### Start by Mao ###
### Oct 3 ###
# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

setwd("C:/PhD/Project/PhD_thesis/mast_trait")

d <- read.csv("data/silvicsClean.csv")

str(d)
summary(d)


# Check missing values
colSums(is.na(d))

d[d == ""] <- NA

ave_freq <- function(x) {
  x <- gsub(" to ", " to ", x) 
  if (grepl(" to ", x)) {
    parts <- as.numeric(strsplit(x, " to ")[[1]])
    return(mean(parts))
  } else {
    return(as.numeric(x))
  }
}

# Make a new column for mast cycle
d$mastCycleAve <- sapply(as.character(d$mastCycle), ave_freq)
head(d$mastCycleAve)

# Make a new column for ave fruit size
d$fruitSizeAve <- sapply(as.character(d$fruitSize.cm.), ave_freq)
d$seedSizeAve <- sapply(as.character(d$seedSize.mm.), ave_freq)

d$latbi <- paste(d$genusName, d$speciesName, sep = "_")

# Monoecious and Dioecious
d$typeMonoOrDio[which(d$typeMonoOrDio == "polygamo-dioecious")] <-"Polygamous"
d$typeMonoOrDio[which(d$typeMonoOrDio == "polygamo-monoecious")] <-"Polygamous"
d$typeMonoOrDio[which(d$typeMonoOrDio == "Polygamo-monoecious")] <-"Polygamous"
d$typeMonoOrDio[which(d$typeMonoOrDio == "Polygamo-dioecious")] <-"Polygamous"
d$typeMonoOrDio[which(d$typeMonoOrDio == "Heterdichogamous")] <-"Polygamous"
d$typeMonoOrDio[which(d$typeMonoOrDio == "Monoecious/dioecious")] <-"Polygamous"


# Masting and pollination
d$pollination[which(d$pollination == "insects")] <- "animals"
d$pollination[which(d$pollination == "wind and insects")] <- "wind and animals"
d$pollination[which(d$pollination == "bird and insects")] <- "animals"
d$pollination[which(d$pollination == "insect and wind")] <- "wind and animals"
d$pollination[which(d$pollination == "insect")] <- "animals"

# Seed dormancy data
d$seedDormancy[is.na(d$seedDormancy)] <- d$dormancyClass[is.na(d$seedDormancy)]
d$seedDormancy[which(d$seedDormancy == "PD")] <- "Y"
d$seedDormancy[which(d$seedDormancy == "ND")] <- "N"
d$seedDormancy[which(d$seedDormancy == "MPD")] <- "Y"
d$seedDormancy[which(d$seedDormancy == "PY")] <- "Y"
d$seedDormancy[which(d$seedDormancy == "PYPD")] <- "Y"
d$seedDormancy[which(d$seedDormancy == "MD")] <- "Y"

# Pollination data
d$pollination[which(d$pollination == "wind, insects")] <- "wind and animals"
d$pollination[which(d$pollination == "bird, insects")] <- "animals"
d$pollination[which(d$pollination == "insects and wind")] <- "wind and animals"

# dispersal data
unique(d$seedDispersal)
d$seedDispersalDetails <- d$seedDispersal
d$seedDispersal[d$seedDispersal %in% c("wind", "wind, water", "fire", "gravity", "gravity, wind", "water")] <- "abiotic"
d$seedDispersalDetails[which(d$seedDispersal == "wind,mammals")] <- "wind, mammals"
d$seedDispersal[d$seedDispersal %in% c("wind,mammals","wind, mammals","gravity, water, birds, mammals","wind, birds","wind, rodents, birds","gravity, wind, rodents, birds","gravity, water, mammals","birds, wind", "water, mammals", "water, birds","gravity, wind, water, birds", "wind, water, birds", "gravity, wind, rodents", "gravity, water, rodents", "water, rodents, birds", "water, mammals, birds", "gravity, wind, mammals", "wind, water, mammals", "gravity, wind, water, mammals")] <- "both"
d$seedDispersal[d$seedDispersal %in% c("mammals", "birds", "birds, mammals", "gravity, birds, mammals", "rodents, birds", "gravity, rodents", "gravity, rodents, birds", "bird", "animals", "gravity, mammals", "mammals, birds", "bird, rodents, mammals", "rodents, mammals", "animals, rodents")] <- "biotic"
d$seedDispersalDetails[which(d$seedDispersal == "bird")] <- "birds"
d$seedDispersalDetails[which(d$seedDispersal == "animals")] <- "mammals"

