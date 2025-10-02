## Started 30 Sep 2025 ##
## By Mao ##

## Get AccSpeciesID for TRY data ##


setwd("C:/PhD/Project/PhD_thesis/mast_trait")
library(dplyr)

d <- read.csv("data/silvics.csv")
# Load TRY data
try_data <- read.delim("data/TryAccSpecies.txt", stringsAsFactors = FALSE)

# Combine genus and species names
d <- d %>%
  mutate(AccSpeciesName = paste(genusName, speciesName, sep = " "))

# Perform the left join to get AccSpeciesID
merged_data <- d %>%
  left_join(try_data, by = "AccSpeciesName")

# Remove missing IDs
matched_ids <- merged_data %>%
  filter(!is.na(AccSpeciesID)) %>%
  pull(AccSpeciesID)

# Collapse into a comma-separated string
id_string <- paste(matched_ids, collapse = ", ")

# Print it
cat(id_string)

# Load downloaded trait data from TRY
try_data <- read.delim("data/TRY/44216.txt", stringsAsFactors = FALSE)

# Subset only leaf longevity data for now
leafLong <- subset(try_data, TraitName %in% "Leaf lifespan (longevity)")
# Order it alphebetically to get a sense what data look like
leafLong <- leafLong[order(leafLong$SpeciesName), ]
# Subset only columns I am interested
leafLong <- leafLong[, c("SpeciesName", "TraitName", "DataName","OrigValueStr","OrigUnitStr","ValueKindName")]
# Get rid of turnover rate data
leafLong <- subset(leafLong, DataName %in% c("Leaf lifespan (longevity, retention time, LL, LLS)","Leaf retention time"))
# Standardize the unit
# Get rid of weird values
leafLong <- subset(leafLong, !(OrigValueStr %in% c(">24", "<=1yr", ">=1yr", "6-12")))
leafLong$OrigValueStr <- as.numeric(leafLong$OrigValueStr)
# Create a new column with values converted to years
leafLong$valueYears <- with(leafLong, ifelse(
  OrigUnitStr %in% c("year", "yr", "years"),
  OrigValueStr,
  ifelse(OrigUnitStr %in% c("month", "months", "mo"),
         OrigValueStr / 12,
         ifelse(OrigUnitStr %in% c("day", "days"),
                OrigValueStr / 365,
                NA))))

# Count how many rows per species
leafLong <- subset(leafLong, lengths(strsplit(SpeciesName, " ")) == 2)
species_counts <- table(leafLong$SpeciesName)

# Identify species with only one row
single_species <- names(species_counts[species_counts == 1])
multi_species <- names(species_counts[species_counts > 1])

# Split dataset
leafLong_single <- leafLong[leafLong$SpeciesName %in% single_species, ]
leafLong_multi <- leafLong[leafLong$SpeciesName %in% multi_species & leafLong$ValueKindName == "Mean", ]

# Combine them
leafLong <- rbind(leafLong_single, leafLong_multi)

# Subset only Species and value
leafLong <- leafLong[, c("SpeciesName","valueYears")]
leafLong <- aggregate(valueYears ~ SpeciesName, data = leafLong, FUN = mean)

names(leafLong) <- c("latbi", "leafLongevity")
d$latbi <- paste(d$genusName, d$speciesName, sep = " ")
d <- d %>%
  left_join(leafLong, by = "latbi")

d$leafLongevity <- ifelse(
  is.na(d$leafLongevity.y),
  d$leafLongevity.x,
  d$leafLongevity.y
)
d$leafLongevity.x <- NULL
d$leafLongevity.y <- NULL

# Subset seed nutrient data
oil <- subset(try_data, TraitName %in% "Seed oil content per seed mass")

# Order it alphebetically to get a sense what data look like
oil <- oil[order(oil$SpeciesName), ]
# Subset only columns I am interested
oil <- oil[, c("SpeciesName", "TraitName", "DataName","OrigValueStr","OrigUnitStr","ValueKindName")]

oil$genus <- str_split_fixed(oil$SpeciesName, ' ', 3)[,1]
oil$species <- str_split_fixed(oil$SpeciesName, ' ', 3)[,2]
oil$latbi <- paste(oil$genus, oil$species, sep = " ")
oil$genus <- NULL
oil$species <- NULL
# Subset only Species and value
oil <- oil[, c("latbi","OrigValueStr")]
oil$OrigValueStr <- as.numeric(oil$OrigValueStr)
oil <- aggregate(OrigValueStr ~ latbi, data = oil, FUN = mean)
names(oil) <- c("latbi", "oilContent")

d <- d %>%
  left_join(oil, by = "latbi")

d$oilContent <- ifelse(
  is.na(d$oilContent.y),
  d$oilContent.x,
  d$oilContent.y
)
d$oilContent.x <- NULL
d$oilContent.y <- NULL

# Subset seed nutrient data
protein <- subset(try_data, TraitName %in% "Seed protein content per seed mass")

# Order it alphebetically to get a sense what data look like
protein <- protein[order(protein$SpeciesName), ]
# Subset only columns I am interested
protein <- protein[, c("SpeciesName", "TraitName", "DataName","OrigValueStr","OrigUnitStr","ValueKindName")]

protein$genus <- str_split_fixed(protein$SpeciesName, ' ', 3)[,1]
protein$species <- str_split_fixed(protein$SpeciesName, ' ', 3)[,2]
protein$latbi <- paste(protein$genus, protein$species, sep = " ")
protein$genus <- NULL
protein$species <- NULL
# Subset only Species and value
protein <- protein[, c("latbi","OrigValueStr")]
protein$OrigValueStr <- as.numeric(protein$OrigValueStr)
protein <- aggregate(OrigValueStr ~ latbi, data = protein, FUN = mean)
names(protein) <- c("latbi", "proteinContent")

d <- d %>%
  left_join(protein, by = "latbi")

d$proteinContent <- ifelse(
  is.na(d$proteinContent.y),
  d$proteinContent.x,
  d$proteinContent.y
)
d$proteinContent.x <- NULL
d$proteinContent.y <- NULL

# Subset seed bank longevity data
seedBank <- subset(try_data, TraitName %in% "Seed (seedbank) longevity")

# Order it alphebetically to get a sense what data look like
seedBank <- seedBank[order(seedBank$SpeciesName), ]
# Subset only columns I am interested
seedBank <- seedBank[, c("SpeciesName", "TraitName", "DataName","OriglName","OrigValueStr","OrigUnitStr")]

d$LongtermSeedBank[which(d$latbi == "Abies grandis")] <-"No"
d$LongtermSeedBank[which(d$latbi == "Abies procera")] <-"No"
d$LongtermSeedBank[which(d$latbi == "Ailanthus altissima")] <-"No"
d$LongtermSeedBank[which(d$latbi == "Alnus glutinosa")] <-"No"
d$LongtermSeedBank[which(d$latbi == "Betula papyrifera")] <-"No"
d$LongtermSeedBank[which(d$latbi == "Betula uber")] <-"No"
d$LongtermSeedBank[which(d$latbi == "Cordia alliodora")] <-"No"
d$LongtermSeedBank[which(d$latbi == "Liriodendron tulipifera")] <-"No"
d$LongtermSeedBank[which(d$latbi == "Melaleuca quinquenervia")] <-"Yes"
d$LongtermSeedBank[which(d$latbi == "Picea glauca")] <-"No"
d$LongtermSeedBank[which(d$latbi == "Picea sitchensis")] <-"No"
d$LongtermSeedBank[which(d$latbi == "Pinus contorta")] <-"Yes"
d$LongtermSeedBank[which(d$latbi == "Pinus monticola")] <-"No"
d$LongtermSeedBank[which(d$latbi == "Pinus nigra")] <-"No"
d$LongtermSeedBank[which(d$latbi == "Pinus strobus")] <-"No"
d$LongtermSeedBank[which(d$latbi == "Pinus sylvestris")] <-"No"
d$LongtermSeedBank[which(d$latbi == "Pinus taeda")] <-"No"
d$LongtermSeedBank[which(d$latbi == "Prunus serotina")] <-"Yes"
d$LongtermSeedBank[which(d$latbi == "Pseudotsuga menziesii")] <-"Yes"
d$LongtermSeedBank[which(d$latbi == "Robinia pseudacacia")] <-"Yes"
d$LongtermSeedBank[which(d$latbi == "Thuja plicata")] <-"Yes"
d$LongtermSeedBank[which(d$latbi == "Tsuga heterophylla")] <-"Yes"

write.csv(d,"silvicsClean.csv")
