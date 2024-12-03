# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

set.seed(42)

n_species <- 10

species <- c(paste("Masting_Species", 1:5), paste("Non_Masting_Species", 1:5))

# Define reproductive traits
# Masting species will have higher average values for traits
flowering_duration_masting <- rnorm(5, mean = 4, sd = 1)
flowering_duration_non_masting <- rnorm(5, mean = 2, sd = 1) 
flowering_duration <- c(flowering_duration_masting, flowering_duration_non_masting)

seed_size_masting <- rnorm(5, mean = 3, sd = 0.5)
seed_size_non_masting <- rnorm(5, mean = 1.5, sd = 0.5)
seed_size <- c(seed_size_masting, seed_size_non_masting)

# Masting species are more likely to have seed dormancy
seed_dormancy_masting <- sample(c(0, 1), 5, replace = TRUE, prob = c(0.3, 0.7))
seed_dormancy_non_masting <- sample(c(0, 1), 5, replace = TRUE, prob = c(0.8, 0.2))
seed_dormancy <- c(seed_dormancy_masting, seed_dormancy_non_masting)

# Masting species are more likely to be monoecious
sex_type_masting <- sample(c("Monoecious", "Dioecious"), 5, replace = TRUE, prob = c(0.7, 0.3))
sex_type_non_masting <- sample(c("Monoecious", "Dioecious"), 5, replace = TRUE, prob = c(0.3, 0.7))
sex_type <- c(sex_type_masting, sex_type_non_masting)

mast_status <- c(rep("Masting", 5), rep("Non-Masting", 5))

data <- data.frame(
  Species = species,
  Mast_Status = mast_status,
  Flowering_Duration = flowering_duration,
  Seed_Size = seed_size,
  Seed_Dormancy = seed_dormancy,
  Sex_Type = sex_type
)

print(data)
# Convert categorical variables into numeric (for PCA)
data$Sex_Type <- as.factor(data$Sex_Type)  # Factor for monoecious/dioecious
data$Mast_Status <- as.factor(data$Mast_Status)  # Factor for masting/non-masting

# Create a new data frame with only numeric variables
data_numeric <- data[, c("Flowering_Duration", "Seed_Size", "Seed_Dormancy")]

data_scaled <- scale(data_numeric)

pca_result <- prcomp(data_scaled, center = TRUE, scale. = TRUE)

summary(pca_result)

# Create a data frame with PCA results
pca_data <- data.frame(
  PC1 = pca_result$x[, 1],  # First principal component
  PC2 = pca_result$x[, 2],  # Second principal component
  Mast_Status = data$Mast_Status,
  Sex_Type = data$Sex_Type
)

library(ggplot2)


# Extract PCA loadings
loadings <- pca_result$rotation[, 1:2]  
loadings <- as.data.frame(loadings)

loadings$Trait <- rownames(loadings)

#loadings <- pca_result$rotation[, 1:2]  # Get the first two principal components loadings
loadings <- as.data.frame(loadings)

# Create a data frame for the traits
loadings$Trait <- rownames(loadings)

ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = Sex_Type, color = Mast_Status),size=4) +  # Plot species points
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1 * 2, yend = PC2 * 2), 
               arrow = arrow(type = "closed", length = unit(0.2, "inches")), 
               size = 1, color = "black") +  # Add arrows indicating trait vectors
  geom_text(data = loadings, aes(x = PC1 * 2.5, y = PC2 * 2.5, label = Trait), 
            size = 3.5, color = "black", vjust = 1.5) +  # Add labels for traits
  labs(title = "PCA Plot with Trait Vectors",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal() +
  scale_color_manual(values = c("Masting" = "blue", "Non-Masting" = "red")) +
  scale_shape_manual(values = c(16, 17)) +
  theme(legend.title = element_blank(), plot.title = element_blank())
