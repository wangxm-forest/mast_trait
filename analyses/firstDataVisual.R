####Data exploration####
####Start by Mao####
####July 13####

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

setwd("C:/PhD/Project/PhD_thesis/mast_trait")

d <- read.csv("data/silvicsClean.csv")

library(ggplot2)
library(corrplot)
library(magrittr)
library(ape)
library(phytools)
library(gridExtra)

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

# Make a subset of conifers only
conifer <- subset(d, familyName %in% c("Pinaceae","Taxodiaceae"))
conifer <- conifer[!is.na(conifer$mastEvent), ]
# Make a subset of angiosperm only
angio <- subset(d, !(familyName %in% c("Pinaceae","Taxodiaceae")))
angio <- angio[!is.na(angio$mastEvent), ]

# Plot seed drop time for masting and non-masting species for conifers
unique(conifer$mainSeedDropTime)


# Mapping from month name to number
month_lookup <- c(
  Jan = 1, Feb = 2, March = 3, April = 4, May = 5, June = 6,
  July = 7, Aug = 8, Sep = 9, Oct = 10, Nov = 11, Dec = 12
)

# Robust parser function
parse_seed_months <- function(time_range) {
  if (is.na(time_range) || time_range == "") return(NA)
  
  # Convert to lowercase and remove extra words
  time_range <- tolower(time_range)
  
  # Extract three-letter month abbreviations
  matches <- regmatches(time_range, gregexpr("\\b[a-z]{3}\\b", time_range))[[1]]
  
  # Capitalize first letter to match lookup names
  matches <- paste0(toupper(substring(matches, 1, 1)), substring(matches, 2, 3))
  
  # Keep only valid months
  valid_months <- matches[matches %in% names(month_lookup)]
  
  if (length(valid_months) < 2) return(NA)  # Need a start and end
  
  # Convert to numeric months
  start_num <- month_lookup[[valid_months[1]]]
  end_num   <- month_lookup[[valid_months[2]]]
  
  if (is.na(start_num) || is.na(end_num)) return(NA)
  
  # Create sequence (handle wrap-around)
  if (start_num > end_num) {
    months_seq <- c(start_num:12, 1:end_num)
  } else {
    months_seq <- start_num:end_num
  }
  
  return(months_seq)
}

expanded_seed_df <- data.frame()

for (i in 1:nrow(conifer)) {
  months <- parse_seed_months(d$mainSeedDropTime[i])
  
  if (!is.na(months[1])) {
    temp <- data.frame(
      species = conifer$latbi[i],
      Masting = conifer$mastEvent[i],
      Month = months
    )
    expanded_seed_df <- rbind(expanded_seed_df, temp)
  }
}

# mastEvent should be a factor
expanded_seed_df$Masting <- factor(expanded_seed_df$Masting, levels = c("Y", "N"))

custom_theme <- theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0.25),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.major = element_line(color = "gray85"),
    panel.grid.minor = element_blank()
  )

ggplot(expanded_seed_df, aes(x = Month, fill = Masting)) +
  geom_histogram(binwidth = 1, boundary = 0.5, color = "black", position = "stack") +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  labs(
    title = "Seed Drop Timing by Masting Status",
    x = "Month",
    y = "Number of Species"
  ) +
  scale_fill_manual(values = c("Y" = "#95B958", "N" = "#F4D166")) +
  custom_theme


# Monoecious and Dioecious
angio$typeMonoOrDio[which(angio$typeMonoOrDio == "polygamo-dioecious")] <-"Polygamous"
angio$typeMonoOrDio[which(angio$typeMonoOrDio == "polygamo-monoecious")] <-"Polygamous"
angio$typeMonoOrDio[which(angio$typeMonoOrDio == "Polygamo-monoecious")] <-"Polygamous"
angio$typeMonoOrDio[which(angio$typeMonoOrDio == "Polygamo-dioecious")] <-"Polygamous"


angio$typeMonoOrDio <- tolower(angio$typeMonoOrDio)
angio$typeMonoOrDio <- ifelse(angio$typeMonoOrDio == "monoecious", "Monoecious",
                                ifelse(angio$typeMonoOrDio == "dioecious", "Dioecious", 
                                       ifelse(angio$typeMonoOrDio == "polygamous", "Polygamous",NA)))
par(mfrow = c(2, 2))
angio$typeMonoOrDio <- factor(angio$typeMonoOrDio, levels = c("Monoecious", "Dioecious","Polygamous"))


mono <- ggplot(angio, aes(x = mastEvent, fill = typeMonoOrDio)) +
  geom_bar(position = "fill", color = "black", width = 0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("Monoecious" = "#95B958", "Dioecious" = "#F4D166","Polygamous" = "lightblue")) +
  labs(
    title = "Reproductive Systems",
    x = "Masting",
    y = "Proportion",
    fill = "Reproductive System"
  ) +
  custom_theme

# Yes, masting species are more likely to be monoecious

# Masting and pollination
angio$pollination[which(angio$pollination == "insects")] <- "animals"
angio$pollination[which(angio$pollination == "wind and insects")] <- "wind and animals"
angio$pollination[which(angio$pollination == "bird and insects")] <- "animals"
angio$pollination[which(angio$pollination == "insect and wind")] <- "wind and animals"
angio$pollination[which(angio$pollination == "insect")] <- "animals"


angio$pollination <- factor(angio$pollination, levels = c("wind", "animals", "wind and animals"))


pollination <- ggplot(angio, aes(x = mastEvent, fill = pollination)) +
  geom_bar(position = "fill", color = "black", width=0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("wind" = "#95B958", "animals" = "#F4D166", "wind and animals" = "#6194BF")) +
  labs(
    title = "Pollination",
    x = "Masting",
    y = "Number of Species",
    fill = "Pollination"  ) +
  custom_theme

# Yes, masting species are more likely to be wind pollinated

# Then let's check on the dormancy and masting
unique(angio$dormancyClass)
angio$dormancyClass[which(angio$dormancyClass == "PD")] <- "Y"
angio$dormancyClass[which(angio$dormancyClass == "ND")] <- "N"
angio$dormancyClass[which(angio$dormancyClass == "MPD")] <- "Y"
angio$dormancyClass[which(angio$dormancyClass == "PY")] <- "Y"
angio$dormancyClass[which(angio$dormancyClass == "PYPD")] <- "Y"
angio$dormancyClass[which(angio$dormancyClass == "MD")] <- "Y"
angio$dormancyClass <- factor(angio$dormancyClass, levels = c("Y", "N"))


dormancy <- ggplot(angio, aes(x = mastEvent, fill = dormancyClass)) +
  geom_bar(position = "fill", color = "black", width = 0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("Y" = "#95B958", "N" = "#F4D166")) +
  labs(
    title = "Seed Dormancy",
    x = "Masting",
    y = "Number of Species",
    fill = "Seed Dormancy" ) +
  custom_theme

# Yes, masting species are more likely to have dormant seed? but there are so many NAs for non-masting species

# Then let's check on the shade tolerance and masting
angio$droughtTolerance <- factor(angio$droughtTolerance, levels = c("High", "Moderate","Low"))

drought <- ggplot(angio, aes(x = mastEvent, fill = droughtTolerance)) +
  geom_bar(position = "fill", color = "black", width = 0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("High" = "#95B958", "Low" = "#F4D166", "Moderate" = "#6194BF")) +
  labs(
    title = "Drought Tolerance",
    x = "Masting",
    y = "Number of Species",
    fill = "Drought Tolerance" ) +
  custom_theme

grid.arrange(mono, pollination, dormancy, drought, nrow = 2, ncol = 2)
# No big difference

# Seed weight


angio$logSeedWeights <- log10(angio$seedWeights)

weight <- ggplot(angio, aes(x = logSeedWeights, fill = mastEvent)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Y" = "#95B958", "N" = "#F4D166")) +
  labs(
    x = "Seed Weight (log10)",
    y = "Density"
  ) +
  custom_theme

angio$logSeedSize <- log10(angio$seedSizeAve)

seedsize <- ggplot(angio, aes(x = logSeedSize, fill = mastEvent)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Y" = "#95B958", "N" = "#F4D166")) +
  labs(
    x = "Seed Size (log10)",
    y = "Density"
  ) +
  custom_theme

angio$logFruitSize <- log10(angio$fruitSizeAve)

fruitsize <- ggplot(angio, aes(x = logFruitSize, fill = mastEvent)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Y" = "#95B958", "N" = "#F4D166")) +
  labs(
    x = "Seed Size (log10)",
    y = "Density"
  ) +
  custom_theme

leafLongevity <- ggplot(angio, aes(x = leafLongevity, fill = mastEvent)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Y" = "#95B958", "N" = "#F4D166")) +
  labs(
    x = "Leaf Longevity",
    y = "Density"
  ) +
  custom_theme

grid.arrange(weight, seedsize, fruitsize, leafLongevity, nrow = 2, ncol = 2)

# Seed weight and mastFrequency
weight <- ggplot(angio, aes(x = logSeedWeights, y = mastCycleAve)) +
  geom_point(alpha = 0.6, color = "#95B958", size = 2) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  labs(
    x = "Seed Weight (log10 scale)",
    y = "Average Mast Frequency (years)"
  ) +
  custom_theme

seedsize <- ggplot(angio, aes(x = logSeedSize, y = mastCycleAve)) +
  geom_point(alpha = 0.6, color = "#95B958", size = 2) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  labs(
    x = "Seed Size (log10 scale)",
    y = "Average Mast Frequency (years)"
  ) +
  custom_theme

fruitsize <- ggplot(angio, aes(x = logFruitSize, y = mastCycleAve)) +
  geom_point(alpha = 0.6, color = "#95B958", size = 2) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  labs(
    x = "Fruit Size (log10 scale)",
    y = "Average Mast Frequency (years)"
  ) +
  custom_theme

leafLongevity <- ggplot(angio, aes(x = leafLongevity, y = mastCycleAve)) +
  geom_point(alpha = 0.6, color = "#95B958", size = 2) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  labs(
    x = "Leaf Longevity",
    y = "Average Mast Frequency (years)"
  ) +
  custom_theme

grid.arrange(weight, seedsize, fruitsize, leafLongevity, nrow = 2, ncol = 2)


# To be honest, I don't really see a pattern here...

#### Conifers ####

# Monoecious and Dioecious
conifer$typeMonoOrDio[which(conifer$typeMonoOrDio == "polygamo-dioecious")] <-"Polygamous"
conifer$typeMonoOrDio[which(conifer$typeMonoOrDio == "polygamo-monoecious")] <-"Polygamous"
conifer$typeMonoOrDio[which(conifer$typeMonoOrDio == "Polygamo-monoecious")] <-"Polygamous"
conifer$typeMonoOrDio[which(conifer$typeMonoOrDio == "Polygamo-dioecious")] <-"Polygamous"


conifer$typeMonoOrDio <- tolower(conifer$typeMonoOrDio)
conifer$typeMonoOrDio <- ifelse(conifer$typeMonoOrDio == "monoecious", "Monoecious",
                              ifelse(conifer$typeMonoOrDio == "dioecious", "Dioecious", 
                                     ifelse(conifer$typeMonoOrDio == "polygamous", "Polygamous",NA)))
par(mfrow = c(2, 2))
conifer$typeMonoOrDio <- factor(conifer$typeMonoOrDio, levels = c("Monoecious", "Dioecious","Polygamous"))


mono <- ggplot(conifer, aes(x = mastEvent, fill = typeMonoOrDio)) +
  geom_bar(position = "fill", color = "black", width = 0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("Monoecious" = "#8EB3D2", "Dioecious" = "#FFB2B6","Polygamous" = "lightblue")) +
  labs(
    title = "Reproductive Systems",
    x = "Masting",
    y = "Proportion",
    fill = "Reproductive System"
  ) +
  custom_theme

# Yes, masting species are more likely to be monoecious

# Masting and pollination
conifer$pollination[which(conifer$pollination == "insects")] <- "animals"
conifer$pollination[which(conifer$pollination == "wind and insects")] <- "wind and animals"
conifer$pollination[which(conifer$pollination == "bird and insects")] <- "animals"
conifer$pollination[which(conifer$pollination == "insect and wind")] <- "wind and animals"
conifer$pollination[which(conifer$pollination == "insect")] <- "animals"


conifer$pollination <- factor(conifer$pollination, levels = c("wind", "animals", "wind and animals"))


pollination <- ggplot(conifer, aes(x = mastEvent, fill = pollination)) +
  geom_bar(position = "fill", color = "black", width=0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("wind" = "#8EB3D2", "animals" = "#FFB2B6", "wind and animals" = "#F4D166")) +
  labs(
    title = "Pollination",
    x = "Masting",
    y = "Number of Species",
    fill = "Pollination"  ) +
  custom_theme

# Yes, masting species are more likely to be wind pollinated

# Then let's check on the dormancy and masting
unique(conifer$dormancyClass)
conifer$dormancyClass[which(conifer$dormancyClass == "PD")] <- "Y"
conifer$dormancyClass[which(conifer$dormancyClass == "ND")] <- "N"
conifer$dormancyClass[which(conifer$dormancyClass == "MPD")] <- "Y"
conifer$dormancyClass[which(conifer$dormancyClass == "PY")] <- "Y"
conifer$dormancyClass[which(conifer$dormancyClass == "PYPD")] <- "Y"
conifer$dormancyClass[which(conifer$dormancyClass == "MD")] <- "Y"
conifer$dormancyClass <- factor(conifer$dormancyClass, levels = c("Y", "N"))


dormancy <- ggplot(conifer, aes(x = mastEvent, fill = dormancyClass)) +
  geom_bar(position = "fill", color = "black", width = 0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("Y" = "#8EB3D2", "N" = "#FFB2B6")) +
  labs(
    title = "Seed Dormancy",
    x = "Masting",
    y = "Number of Species",
    fill = "Seed Dormancy" ) +
  custom_theme

# Yes, masting species are more likely to have dormant seed? but there are so many NAs for non-masting species

# Then let's check on the shade tolerance and masting
conifer$droughtTolerance <- factor(conifer$droughtTolerance, levels = c("High", "Moderate","Low"))

drought <- ggplot(conifer, aes(x = mastEvent, fill = droughtTolerance)) +
  geom_bar(position = "fill", color = "black", width = 0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("High" = "#8EB3D2", "Low" = "#FFB2B6", "Moderate" = "#F4D166")) +
  labs(
    title = "Drought Tolerance",
    x = "Masting",
    y = "Number of Species",
    fill = "Drought Tolerance" ) +
  custom_theme

grid.arrange(mono, pollination, dormancy, drought, nrow = 2, ncol = 2)
# No big difference

# Seed weight


conifer$logSeedWeights <- log10(conifer$seedWeights)

weight <- ggplot(conifer, aes(x = logSeedWeights, fill = mastEvent)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Y" = "#8EB3D2", "N" = "#FFB2B6")) +
  labs(
    x = "Seed Weight (log10)",
    y = "Density"
  ) +
  custom_theme

conifer$logSeedSize <- log10(conifer$seedSizeAve)

seedsize <- ggplot(conifer, aes(x = logSeedSize, fill = mastEvent)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Y" = "#8EB3D2", "N" = "#FFB2B6")) +
  labs(
    x = "Seed Size (log10)",
    y = "Density"
  ) +
  custom_theme

conifer$logFruitSize <- log10(conifer$fruitSizeAve)

fruitsize <- ggplot(conifer, aes(x = logFruitSize, fill = mastEvent)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Y" = "#8EB3D2", "N" = "#FFB2B6")) +
  labs(
    x = "Seed Size (log10)",
    y = "Density"
  ) +
  custom_theme

leafLongevity <- ggplot(conifer, aes(x = leafLongevity, fill = mastEvent)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Y" = "#8EB3D2", "N" = "#FFB2B6")) +
  labs(
    x = "Leaf Longevity",
    y = "Density"
  ) +
  custom_theme

grid.arrange(weight, seedsize, fruitsize, leafLongevity, nrow = 2, ncol = 2)

# Seed weight and mastFrequency
weight <- ggplot(conifer, aes(x = logSeedWeights, y = mastCycleAve)) +
  geom_point(alpha = 0.6, color = "#8EB3D2", size = 2) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  labs(
    x = "Seed Weight (log10 scale)",
    y = "Average Mast Frequency (years)"
  ) +
  custom_theme

seedsize <- ggplot(conifer, aes(x = logSeedSize, y = mastCycleAve)) +
  geom_point(alpha = 0.6, color = "#8EB3D2", size = 2) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  labs(
    x = "Seed Size (log10 scale)",
    y = "Average Mast Frequency (years)"
  ) +
  custom_theme

fruitsize <- ggplot(conifer, aes(x = logFruitSize, y = mastCycleAve)) +
  geom_point(alpha = 0.6, color = "#8EB3D2", size = 2) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  labs(
    x = "Fruit Size (log10 scale)",
    y = "Average Mast Frequency (years)"
  ) +
  custom_theme

leafLongevity <- ggplot(conifer, aes(x = leafLongevity, y = mastCycleAve)) +
  geom_point(alpha = 0.6, color = "#8EB3D2", size = 2) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  labs(
    x = "Leaf Longevity",
    y = "Average Mast Frequency (years)"
  ) +
  custom_theme

grid.arrange(weight, seedsize, fruitsize, leafLongevity, nrow = 2, ncol = 2)
