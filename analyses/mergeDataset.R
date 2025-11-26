#################################################
#########Started by Mao on Nov 25, 2025##########
#################################################

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

# working directory
setwd("C:/PhD/Project/PhD_thesis/mast_trait")

# read in MASTREE data
mastree <- read.csv("C:/PhD/Project/MASTREEplus/Data/MASTREEplus_2024-06-26_V2.csv",header = TRUE)
mastree <- subset(mastree, 
                      select = c(Species,VarType,Variable, Value,Start,End))

# read in masting change over time data
masting <- read.csv("C:/PhD/Project/MASTREEplus/Data/masting_change_over_time.csv", header = TRUE)
masting <- subset(masting, 
                  select = c(SPECIES,m_CV))
masting$SPECIES <- sub("_", " ", masting$SPECIES)
masting$SPECIES <- str_to_sentence(masting$SPECIES)

# read in silvics data
silvics <- read.csv("data/cleanSilvics.csv", header = TRUE)

silvics$latbi[!is.na(silvics$updateName)] <- silvics$updateName[!is.na(silvics$updateName)]
silvics$latbi <- sub("_", " ", silvics$latbi)


silvics_overlap <- mastree[mastree$Species %in% silvics$latbi, ]
silvics_overlap1 <- masting[masting$SPECIES %in% silvics$latbi, ]
silvics_overlap1 <- distinct(silvics_overlap1)
names(silvics_overlap1)[1]<-paste("latbi")

write.csv(silvics_overlap,"output/silvicsMASTREE.csv")
silvics <- left_join(silvics, silvics_overlap1, by="latbi")
write.csv(silvics, "data/cleanSilvics.csv")
