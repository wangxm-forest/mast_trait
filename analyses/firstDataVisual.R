####Data exploration####
####Start by Mao####
####July 13####

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

setwd("C:/PhD/Project/PhD_thesis/mast_trait")

d <- read.csv("data/silvics1.csv")

library(ggplot2)
library(corrplot)
library(rpart)
library(rpart.plot)
library(magrittr)

str(d)
summary(d)


# Check missing values
colSums(is.na(d))


# Combine genus and species
d$latbi <- paste(d$genusName, d$speciesName, sep = "_")

# Make a new column for mast cycle
ave_freq <- function(x) {
  x <- gsub(" to ", " to ", x)  # Replace en-dash if needed
  if (grepl(" to ", x)) {
    parts <- as.numeric(strsplit(x, " to ")[[1]])
    return(mean(parts))
  } else {
    return(as.numeric(x))
  }
}

d$mastCycleAve <- sapply(as.character(d$mastCycle), ave_freq)
head(d$mastCycleAve)

# Plot seed drop time for masting and non-masting species

unique(d$mainSeedDropTime)


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

for (i in 1:nrow(d)) {
  months <- parse_seed_months(d$mainSeedDropTime[i])
  
  if (!is.na(months[1])) {
    temp <- data.frame(
      species = d$latbi[i],
      Masting = d$mastEvent[i],
      Month = months
    )
    expanded_seed_df <- rbind(expanded_seed_df, temp)
  }
}

# mastEvent should be a factor
expanded_seed_df$Masting <- factor(expanded_seed_df$Masting, levels = c("Y", "N"))

custom_theme <- theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
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
d$typeMonoOrDio <- tolower(d$typeMonoOrDio)
d$typeMonoOrDio <- ifelse(d$typeMonoOrDio == "monoecious", "Monoecious",
                                ifelse(d$typeMonoOrDio == "dioecious", "Dioecious", NA))
d$typeMonoOrDio <- factor(d$typeMonoOrDio, levels = c("Monoecious", "Dioecious"))

ggplot(d, aes(x = mastEvent, fill = typeMonoOrDio)) +
  geom_bar(color = "black") +
  scale_fill_manual(values = c("Monoecious" = "#95B958", "Dioecious" = "#F4D166")) +
  labs(
    title = "Reproductive System of Masting vs Non-Masting Species",
    x = "Masting",
    y = "Number of Species",
    fill = "Reproductive System"
  ) +
  custom_theme


ggplot(d[!is.na(d$mastEvent) & d$mastEvent != "", ], aes(x = mastEvent, fill = typeMonoOrDio)) +
  geom_bar(position = "fill", color = "black", width = 0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("Monoecious" = "#95B958", "Dioecious" = "#F4D166")) +
  labs(
    title = "Proportion of Reproductive Systems in Masting vs Non-Masting Species",
    x = "Masting",
    y = "Proportion",
    fill = "Reproductive System"
  ) +
  custom_theme

# Yes, masting species are more likely to be monoecious

# Masting and pollination
d$pollination[which(d$pollination == "insects")] <- "animals"
d$pollination[which(d$pollination == "wind and insects")] <- "wind and animals"
d$pollination[which(d$pollination == "bird and insects")] <- "animals"
d$pollination[which(d$pollination == "insect and wind")] <- "wind and animals"
d$pollination[which(d$pollination == "insect")] <- "animals"


d$pollination <- factor(d$pollination, levels = c("wind", "animals", "wind and animals"))

ggplot(d[!is.na(d$mastEvent) & d$mastEvent != "", ], aes(x = mastEvent, fill = pollination)) +
  geom_bar(color = "black", width = 0.4) +
  scale_fill_manual(values = c("wind" = "#95B958", "animals" = "#F4D166", "wind and animals" = "#6194BF")) +
  labs(
    title = "Pollination of Masting vs Non-Masting Species",
    x = "Masting",
    y = "Number of Species",
    fill = "Pollination"
  ) +
  custom_theme

ggplot(d[!is.na(d$mastEvent) & d$mastEvent != "", ], aes(x = mastEvent, fill = pollination)) +
  geom_bar(position = "fill", color = "black", width=0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("wind" = "#95B958", "animals" = "#F4D166", "wind and animals" = "#6194BF")) +
  labs(
    title = "Pollination of Masting vs Non-Masting Species",
    x = "Masting",
    y = "Number of Species",
    fill = "Pollination"  ) +
  custom_theme

# Yes, masting species are more likely to be wind pollinated

# Then let's check on the dormancy and masting
unique(d$dormancyClass)
d$dormancyClass[which(d$dormancyClass == "PD")] <- "Y"
d$dormancyClass[which(d$dormancyClass == "ND")] <- "N"
d$dormancyClass[which(d$dormancyClass == "MPD")] <- "Y"
d$dormancyClass[which(d$dormancyClass == "PY")] <- "Y"
d$dormancyClass[which(d$dormancyClass == "PYPD")] <- "Y"
d$dormancyClass[which(d$dormancyClass == "MD")] <- "Y"
d$dormancyClass <- factor(d$dormancyClass, levels = c("Y", "N"))


ggplot(d, aes(x = mastEvent, fill = dormancyClass)) +
  geom_bar(position = "fill", color = "black", width = 0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("Y" = "#95B958", "N" = "#F4D166")) +
  labs(
    title = "Seed Dormancy of Masting vs Non-Masting Species",
    x = "Masting",
    y = "Number of Species",
    fill = "Seed Dormancy" ) +
  custom_theme

# Yes, masting species are more likely to have dormant seed? but there are so many NAs for non-masting species

# Then let's check on the dormancy and masting
d$shadeTolerance <- factor(d$shadeTolerance, levels = c("tolerant", "intolerant","intermediate"))

ggplot(d, aes(x = mastEvent, fill = shadeTolerance)) +
  geom_bar(color = "black") +
  scale_fill_manual(values = c("tolerant" = "#95B958", "intolerant" = "#F4D166", "intermediate" = "#6194BF")) +
  labs(
    title = "Shade Tolerance of Masting vs Non-Masting Species",
    x = "Masting",
    y = "Number of Species",
    fill = "Shade Tolerance"
  ) +
  theme_minimal()

ggplot(d[!is.na(d$mastEvent) & d$mastEvent != "", ], aes(x = mastEvent, fill = shadeTolerance)) +
  geom_bar(position = "fill", color = "black", width = 0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("tolerant" = "#95B958", "intolerant" = "#F4D166", "intermediate" = "#6194BF")) +
  labs(
    title = "Shade Tolerance of Masting vs Non-Masting Species",
    x = "Masting",
    y = "Number of Species",
    fill = "Shade Tolerance" ) +
  custom_theme

# No big difference

# Seed weight


d$logSeedWeights <- log10(d$seedWeights)

ggplot(d, aes(x = logSeedWeights, fill = mastEvent)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Y" = "#95B958", "N" = "#F4D166")) +
  labs(
    title = "Density of Seed Weight (Log Scale)",
    x = "Seed Weight (log10)",
    y = "Density"
  ) +
  custom_theme

# Seed weight and mastFrequency


ggplot(d, aes(x = logSeedWeights, y = mastCycleAve)) +
  geom_point(alpha = 0.6, color = "#95B958", size = 2) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  labs(
    title = "Seed Weight vs Average Mast Frequency (Log Scale)",
    x = "Seed Weight (log10 scale)",
    y = "Average Mast Frequency (years)"
  ) +
  custom_theme

# To be honest, I don't really see a pattern here...


