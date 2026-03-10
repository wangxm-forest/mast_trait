#### Phylogeny logistic regression ####
## Started by Mao ##
## Nov-28-2025 ##

library(ape)
library(caper)
library(geiger)
library(phylolm)
library(brglm2)
library(detectseparation)
library(gridExtra)
library(grid)
library(dplyr)
library(ggplot2)
library(xtable)
library(patchwork)

rm(list = ls())
options(stringsAsFactors = FALSE)

setwd("C:/PhD/Project/PhD_thesis/mast_trait")

## Data preparation ----

# Load & prepare Data
d <- read.csv("data/cleanSilvics.csv")
phytree <- read.tree("output/silvicsPhylogenyFull.tre")

d$latbi <- gsub(" ", "_", d$latbi)
d <- d[!is.na(d$mastEvent), ]
d$mastEvent <- ifelse(d$mastEvent == "Y", 1, 0)

# Create group variable
d$group <- ifelse(d$familyName %in% c("Pinaceae", "Taxodiaceae"),
                  "conifer", "angiosperm")
d$group <- factor(d$group)

conifer <- d[d$familyName %in% c("Pinaceae", "Taxodiaceae"), ]
angio   <- d[!(d$familyName %in% c("Pinaceae", "Taxodiaceae")), ]


# Factorize all the categorical traits, I am doing this separately because conifer and angio have different levels for some traits
conifer$droughtTolerance <- as.factor(conifer$droughtTolerance)
conifer$typeMonoOrDio <- as.factor(conifer$typeMonoOrDio)
conifer$pollination <- as.factor(conifer$pollination)
conifer$seedDispersal <- as.factor(conifer$seedDispersal)
conifer$seedDormancy <- as.factor(conifer$seedDormancy)
angio$droughtTolerance <- as.factor(angio$droughtTolerance)
angio$typeMonoOrDio <- as.factor(angio$typeMonoOrDio)
angio$pollination <- as.factor(angio$pollination)
angio$seedDispersal <- as.factor(angio$seedDispersal)
angio$seedDormancy <- as.factor(angio$seedDormancy)

# log10-transform + scale
#log_scale <- function(x) scale(log10(x))

# Continuous traits
d$logSeedWeight <- log(d$seedWeights)
conifer$logSeedWeight <- log(conifer$seedWeights)
angio$logSeedWeight   <- log(angio$seedWeights)

d$logFruit <- log(d$fruitSizeAve)
conifer$logFruit <- log(conifer$fruitSizeAve)
angio$logFruit <- log(angio$fruitSizeAve)

d$logSeedSize <- log(d$seedSizeAve)
conifer$logSeedSize <- log(conifer$seedSizeAve)
angio$logSeedSize   <- log(angio$seedSizeAve)

phytree$node.label <- NULL

phyconifer <- drop.tip(phytree, setdiff(phytree$tip.label, conifer$latbi))
phyangio   <- drop.tip(phytree, setdiff(phytree$tip.label, angio$latbi))

rownames(conifer) <- conifer$latbi
rownames(angio)   <- angio$latbi



## Functions used later for the analysis ----

# Extract key results from a phyloglm object
tidy_phyloglm <- function(model) {
  s <- summary(model)
  coefs <- s$coefficients
  N <- nrow(model$X)
  out <- data.frame(
    term = rownames(coefs),
    estimate = coefs[, 1],
    std_error = coefs[, 2],
    z_value = coefs[, 3],
    p_value = coefs[, 4],
    alpha = model$alpha,
    N = N,
    row.names = NULL
  )
  return(out)
}

# Run phyloglm model and attach metadata
run_phyloglm <- function(formula, data, phy, method, trait_name) {
  model <- phyloglm(formula, phy = phy, data = data, method = method)
  tbl <- tidy_phyloglm(model)
  tbl$trait <- trait_name
  tbl$model_formula <- deparse(formula)
  return(tbl)
}

# Get legend
get_legend <- function(myplot) {
  g <- ggplotGrob(myplot)
  leg <- g$grobs[which(sapply(g$grobs, function(x) x$name) == "guide-box")]
  leg[[1]]
}

# Model definitions

conifer_list <- list(
  list(name="Seed dispersal", formula=mastEvent ~ seedDispersal, data=conifer, phy=phyconifer, method="logistic_MPLE"),
  list(name="Seed dormancy",  formula=mastEvent ~ seedDormancy, data=conifer, phy=phyconifer, method="logistic_MPLE"),
  list(name="Reproductive type",       formula=mastEvent ~ typeMonoOrDio, data=conifer, phy=phyconifer, method="logistic_MPLE"),
  list(name="Seed weight (log)",    formula=mastEvent ~ logSeedWeight, data=conifer, phy=phyconifer, method="logistic_IG10"),
  list(name="Fruit size (log)",     formula=mastEvent ~ logFruit, data=conifer, phy=phyconifer, method="logistic_IG10"),
  list(name="Seed size (log)",      formula=mastEvent ~ logSeedSize, data=conifer, phy=phyconifer, method="logistic_IG10"),
  list(name="Leaf longevity (years)", formula=mastEvent ~ leafLongevity, data=conifer, phy=phyconifer, method="logistic_IG10"),
  list(name="Drought tolerance",    formula=mastEvent ~ droughtTolerance, data=conifer, phy=phyconifer, method="logistic_IG10")
)

angio_list <- list(
list(name="Dispersal mode",   formula=mastEvent ~ seedDispersal, data=angio,   phy=phyangio,   method="logistic_MPLE"),
list(name="Pollination mode",   formula=mastEvent ~ pollination, data=angio,   phy=phyangio,   method="logistic_MPLE"),
list(name="Seed dormancy",    formula=mastEvent ~ seedDormancy, data=angio,   phy=phyangio,   method="logistic_MPLE"),
list(name="Reproductive type",         formula=mastEvent ~ typeMonoOrDio, data=angio,   phy=phyangio,   method="logistic_MPLE"),
list(name="Seed weight (log)",      formula=mastEvent ~ logSeedWeight, data=angio,   phy=phyangio,   method="logistic_IG10"),
list(name="Fruit size (log)",       formula=mastEvent ~ logFruit, data=angio,   phy=phyangio,   method="logistic_IG10"),
list(name="Seed size (log)",        formula=mastEvent ~ logSeedSize, data=angio,   phy=phyangio,   method="logistic_IG10"),
list(name="Oil content (%)",      formula=mastEvent ~ oilContent, data=angio, phy=phyangio, method="logistic_IG10"),
list(name="Leaf longevity (years)",   formula=mastEvent ~ leafLongevity, data=angio,   phy=phyangio,   method="logistic_IG10"),
list(name="Drought tolerance",      formula=mastEvent ~ droughtTolerance, data=angio,   phy=phyangio,   method="logistic_IG10")
)


# Run models

conifer_results <- NULL
angio_results <- NULL

for (m in conifer_list) {
  cat("Running model:", m$name, "\n")
  
  tbl <- run_phyloglm(
    formula = m$formula,
    data    = m$data,
    phy     = m$phy,
    method  = m$method,
    trait_name = m$name
  )
  
  conifer_results <- rbind(conifer_results, tbl)
}

for (m in angio_list) {
  cat("Running model:", m$name, "\n")
  
  tbl <- run_phyloglm(
    formula = m$formula,
    data    = m$data,
    phy     = m$phy,
    method  = m$method,
    trait_name = m$name
  )
  
  angio_results <- rbind(angio_results, tbl)
}

#Make a table to present the results:
clean_results <- function(results) {
  
  results_no_int <- results[results$term != "(Intercept)", ]
  
  results_no_int$signif <- cut(
    results_no_int$p_value,
    breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
    labels = c("***", "**", "*", ".", "")
  )
  
  results_no_int$estimate  <- round(results_no_int$estimate,  3)
  results_no_int$std_error <- round(results_no_int$std_error, 3)
  results_no_int$z_value   <- round(results_no_int$z_value,   2)
  results_no_int$p_value   <- signif(results_no_int$p_value,  3)
  results_no_int$alpha     <- round(results_no_int$alpha,     3)
  
  final_table <- results_no_int[, c(
    "trait",
    "term",
    "estimate",
    "std_error",
    "p_value",
    "alpha",
    "N"
  )]
  
  colnames(final_table) <- c(
    "Trait",
    "Predictor",
    "Estimate",
    "SE",
    "P",
    "Phylo α",
    "N"
  )
  final_table$Predictor <- gsub("logFruitStd",     "Fruit size (log, std)", final_table$Predictor)
  final_table$Predictor <- gsub("logSeedWeightStd","Seed weight (log, std)", final_table$Predictor)
  final_table$Predictor <- gsub("logSeedSizeStd",  "Seed size (log, std)", final_table$Predictor)
  final_table$Predictor <- gsub("seedDispersalbiotic",  "Biotic vs abiotic", final_table$Predictor)
  final_table$Predictor <- gsub("seedDispersalboth",  "Both vs abiotic", final_table$Predictor)
  final_table$Predictor <- gsub("pollinationwind",  "Wind vs animal pollination", final_table$Predictor)
  final_table$Predictor <- gsub("pollinationwind and animals",  "Wind+animal vs animal pollination", final_table$Predictor)
  final_table$Predictor <- gsub("seedDormancyY",  "Dormant vs non−dormant", final_table$Predictor)
  final_table$Predictor <- gsub("typeMonoOrDioMonoecious",  "Monoecious vs dioecious", final_table$Predictor)
  final_table$Predictor <- gsub("typeMonoOrDioPolygamous",  "Monoecious vs polygamous", final_table$Predictor)
  final_table$Predictor <- gsub("oilContent",  "Oil content (%)", final_table$Predictor)
  final_table$Predictor <- gsub("leafLongevity",  "Leaf longevity (years)", final_table$Predictor)
  final_table$Predictor <- gsub("droughtToleranceLow",  "Low vs high drought tolerance", final_table$Predictor)
  final_table$Predictor <- gsub("droughtToleranceModerate",  "Moderate vs high drought tolerance", final_table$Predictor) 
  return(final_table)
}

conifer_results <- clean_results(conifer_results)
conifer_results <- cbind(Group = "Gymnosperm", conifer_results)

angio_results <- clean_results(angio_results)
angio_results <- cbind(Group = "Angiosperm", angio_results)

d_results <- rbind(conifer_results, angio_results)

dTable <- xtable(d_results, 
               caption = "Phylogenetic Generalized Linear Model Results", 
               label = "tab:regression")
print(dTable, type = "latex", include.rownames = FALSE)




###Plot the raw data with the phyloglm results ----
## Angiosperm
getAnnotation <- function(trait_name, model_results, subData) {
  
  df <- model_results[model_results$Trait == trait_name, ]
  
  # Always keep N
 
  n_text <- paste0("Nmast = ", sum(subData[,1] == "1" & !is.na(subData[,2])), " Nnon = ",  sum(subData[,1] == "0" & !is.na(subData[,2])))
  
  # Keep only predictors with P < 0.5
  df_sig <- df[df$P < 0.05, ]
  alpha_text <- paste0("Phylo α = ", round(unique(df$`Phylo α`), 2))
  # If none meet threshold → return only alpha
  if(nrow(df_sig) == 0){
    return(paste(n_text, sep = "\n"))
  }
  
  # If multiple predictors (categorical trait)
  if(length(unique(df$Predictor)) > 1){
    
    pred_text <- paste(
      df_sig$Predictor,
      ": Estimate = ", round(df_sig$Estimate, 2),
      " ± ", round(df_sig$SE, 2),
      " P=", round(df_sig$P, 2),
      " *",
      sep = "",
      collapse = "\n"
    )
    
    label <- paste(alpha_text, pred_text, n_text, sep = "\n")
    
  } else {
    
    # Single predictor
    pred_text <- paste0(
      " Estimate=", round(df_sig$Estimate, 2),
      " ± ", round(df_sig$SE, 2),
      " P=", round(df_sig$P, 2),
      " *"
    )
    
    label <- paste(alpha_text, pred_text, n_text, sep = "\n")
  }
  
  return(label)
}

getAnnotation("Dispersal mode", angio_results, dispData)
# "n=50, P=0.03,0.12"

dispData <- angio[, c("mastEvent", "seedDispersal")]
dispData$mastEvent <- as.factor(dispData$mastEvent)
dispData$seedDispersal <- as.factor(dispData$seedDispersal)
# Create ggplot object
disp <- ggplot(dispData, aes(x = seedDispersal, fill = mastEvent)) +
  geom_bar(position = "dodge", colour = "black") +
annotate("text", x = 2.5, y = 29, 
         label = getAnnotation("Dispersal mode", angio_results, dispData),
         size = 2)+
  labs(y = "Count", x = NULL) +
  theme_classic(base_size = 12) + scale_x_discrete(labels = c("abiotic" = "Abiotic", "biotic" = "Biotic", "both"="Both")) + scale_fill_manual(
    name = "Masting", 
    values = c("0"="#A9D5B1", "1"="#ED562C"), 
    labels = c("0"="No", "1"="Yes") 
  ) + scale_y_continuous(limits = c(0, 30))
legend <- get_legend(disp)
disp <- ggplot(dispData, aes(x = seedDispersal, fill = mastEvent)) +
  geom_bar(position = "dodge", colour = "black") +
  annotate("text", x = 2.5, y = 29, 
           label = getAnnotation("Dispersal mode", angio_results, dispData),
           size = 2)+
  labs(y = "Count (Seed dispersal)", x = NULL) +
  theme_classic(base_size = 12) + scale_x_discrete(labels = c("abiotic" = "Abiotic", "biotic" = "Biotic", "both"="Both")) + scale_fill_manual(
    name = "Masting", 
    values = c("0"="#A9D5B1", "1"="#ED562C"), 
    labels = c("0"="No", "1"="Yes") 
  ) + scale_y_continuous(limits = c(0, 30)) + theme(legend.position = "none")

dormData <- angio[, c("mastEvent", "seedDormancy")]
dormData$mastEvent <- as.factor(dormData$mastEvent)
dormData$seedDormancy <- as.factor(dormData$seedDormancy)
# Create ggplot object
dorm <- ggplot(dormData, aes(x = seedDormancy, fill = mastEvent)) +
  geom_bar(position = "dodge", colour = "black") +
  annotate("text", x = 2, y = 43, 
           label = getAnnotation("Seed dormancy", angio_results, dormData),
           size = 2) +
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#A9D5B1", "1"="#ED562C"), 
    labels = c("0"="No", "1"="Yes") 
  ) + scale_x_discrete(labels = c("N" = "No dormancy", "Y" = "Dormancy")) +
  labs(y = "Count (Seed dormancy)", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(0, 45)) + theme(legend.position = "none")

weightData <- angio[, c("mastEvent", "logSeedWeight")]
weightData$mastEvent <- as.factor(weightData$mastEvent)

weight <- ggplot(weightData, aes(x = mastEvent, y = logSeedWeight, fill = mastEvent)) +
  geom_violin(trim = FALSE, alpha = 0.5, colour = "black") +
  geom_point(position = position_jitter(width = 0.08), size = 1.5, alpha = 0.7)  +
annotate("text", x = 1.5, y = 17, 
         label = getAnnotation("Seed weight (log)", angio_results, weightData),
         size = 2)+
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#A9D5B1", "1"="#ED562C"), 
    labels = c("0"="No", "1"="Yes")) + scale_x_discrete(labels = c("0" = "Non-masting", "1" = "Masting")) +
  labs(y = "Seed Weight (log)", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(NA, 18))+ theme(legend.position = "none")

fruitData <- angio[, c("mastEvent", "logFruit")]
fruitData$mastEvent <- as.factor(fruitData$mastEvent)

fruit <- ggplot(fruitData, aes(x = mastEvent, y = logFruit, fill = mastEvent)) +
  geom_violin(trim = FALSE, alpha = 0.5, colour = "black") +
annotate("text", x = 1.5, y = 5.7, 
         label = getAnnotation("Fruit size (log)", angio_results, fruitData),
         size = 2) +
  geom_point(position = position_jitter(width = 0.08), size = 1.5, alpha = 0.7) +
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#A9D5B1", "1"="#ED562C"), 
    labels = c("0"="No", "1"="Yes") 
  ) + scale_x_discrete(labels = c("0" = "Non-masting", "1" = "Masting")) +
  labs(y = "Fruit size (log)", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(NA, 6))+ theme(legend.position = "none")

oilData <- angio[, c("mastEvent", "oilContent")]
oilData$mastEvent <- as.factor(oilData$mastEvent)

oil <- ggplot(oilData, aes(x = mastEvent, y = oilContent, fill = mastEvent)) +
  geom_violin(trim = FALSE, alpha = 0.5, colour = "black") +
annotate("text", x = 1.5, y = 115, 
         label = getAnnotation("Oil content (%)", angio_results, oilData),
         size = 2) +
  geom_point(position = position_jitter(width = 0.08), size = 1.5, alpha = 0.7) +
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#A9D5B1", "1"="#ED562C"), 
    labels = c("0"="No", "1"="Yes") 
  ) + scale_x_discrete(labels = c("0" = "Non-masting", "1" = "Masting")) +
  labs(y = "Oil content %", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(NA, 120))+ theme(legend.position = "none")

pollData <- angio[, c("mastEvent", "pollination")]
pollData$mastEvent <- as.factor(pollData$mastEvent)
pollData$pollination <- as.factor(pollData$pollination)
# Create ggplot object
poll <- ggplot(pollData, aes(x = pollination, fill = mastEvent)) +
  geom_bar(position = "dodge", colour = "black")  +
annotate("text", x = 2.5, y = 33, 
         label = getAnnotation("Pollination mode", angio_results, pollData),
         size = 2) +
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#A9D5B1", "1"="#ED562C"), 
    labels = c("0"="No", "1"="Yes") 
  ) + scale_x_discrete(labels = c("animals" = "Animal", "wind" = "Wind", "wind and animals"="Both")) +
  labs(y = "Count (Pollination mode)", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(0, 35))+ theme(legend.position = "none")

repData <- angio[, c("mastEvent", "typeMonoOrDio")]
repData$mastEvent <- as.factor(repData$mastEvent)
repData$typeMonoOrDio <- as.factor(repData$typeMonoOrDio)
# Create ggplot object
rep <- ggplot(repData, aes(x = typeMonoOrDio, fill = mastEvent)) +
  geom_bar(position = "dodge", colour = "black") + annotate("text", x = 2.5, y = 48, 
         label = getAnnotation("Reproductive type", angio_results, repData),
         size = 2)+
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#A9D5B1", "1"="#ED562C"), 
    labels = c("0"="No", "1"="Yes") 
  ) +
  labs(y = "Count (Reproductive type)", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(0, 50))+ theme(legend.position = "none")

droughtData <- angio[, c("mastEvent", "droughtTolerance")]
droughtData$mastEvent <- as.factor(droughtData$mastEvent)
droughtData$droughtTolerance <- as.factor(droughtData$droughtTolerance)
# Create ggplot object
drought <- ggplot(droughtData, aes(x = droughtTolerance, fill = mastEvent)) +
  geom_bar(position = "dodge", colour = "black") + annotate("text", x = 2.5, y = 24.5, label = getAnnotation("Drought tolerance", angio_results, droughtData),
         size = 2) +
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#A9D5B1", "1"="#ED562C"), 
    labels = c("0"="No", "1"="Yes") 
  ) +
  labs(y = "Count (Drought tolerance)", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(0, 25))+ theme(legend.position = "none")

leafData <- angio[, c("mastEvent", "leafLongevity")]
leafData$mastEvent <- as.factor(leafData$mastEvent)

leaf <- ggplot(leafData, aes(x = mastEvent, y = leafLongevity, fill = mastEvent)) +
  geom_violin(trim = FALSE, alpha = 0.5, colour = "black") +
  geom_point(position = position_jitter(width = 0.08), size = 1.5, alpha = 0.7) + annotate("text", x = 1.5, y = 3.85, label = getAnnotation("Leaf longevity (years)", angio_results, leafData),
         size = 2) +
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#A9D5B1", "1"="#ED562C"), 
    labels = c("0"="No", "1"="Yes"))  + scale_x_discrete(labels = c("0" = "Non-masting", "1" = "Masting")) +
  labs(y = "Leaf longevity (year)", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(NA, 4))+ theme(legend.position = "none")

predation <- weight + fruit + oil + disp + dorm + guides(color = "none")
pollination <- poll + rep + plot_layout(ncol = 3, widths = c(1,1,1)) 
resource <- drought + leaf + plot_layout(ncol = 3, widths = c(1,1,1)) + guides(color = "none")

final <- predation / pollination / resource +
  plot_layout(heights = c(1, 0.5, 0.5), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(face = "bold", size = 10)
  )
ggsave(
  filename = "output/figures/rawDataWithStatsAngio.pdf",
  plot = final,
  width = 12,
  height = 12
)

## Gymnosperm
dispData <- conifer[, c("mastEvent", "seedDispersal")]
dispData$mastEvent <- as.factor(dispData$mastEvent)
dispData$seedDispersal <- as.factor(dispData$seedDispersal)
# Create ggplot object
disp <- ggplot(dispData, aes(x = seedDispersal, fill = mastEvent)) +
  geom_bar(position = "dodge", colour = "black") +
annotate("text", x = 2.5, y = 38, label = getAnnotation("Seed dispersal", conifer_results, dispData),
         size = 2) +
  labs(y = "Count (Seed dispersal)", x = NULL) +
  theme_classic(base_size = 12) + scale_x_discrete(labels = c("abiotic" = "Abiotic", "biotic" = "Biotic", "both"="Both")) + scale_fill_manual(
    name = "Masting", 
    values = c("0"="#A9D5B1", "1"="#ED562C"), 
    labels = c("0"="No", "1"="Yes") 
  ) + scale_y_continuous(limits = c(0, 40))

dormData <- conifer[, c("mastEvent", "seedDormancy")]
dormData$mastEvent <- as.factor(dormData$mastEvent)
dormData$seedDormancy <- as.factor(dormData$seedDormancy)
# Create ggplot object
dorm <- ggplot(dormData, aes(x = seedDormancy, fill = mastEvent)) +
  geom_bar(position = "dodge", colour = "black") +
  annotate("text", x = 2, y = 47, label = getAnnotation("Seed dormancy", conifer_results, dormData),
           size = 2) +
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#A9D5B1", "1"="#ED562C"), 
    labels = c("0"="No", "1"="Yes") 
  ) + scale_x_discrete(labels = c("N" = "No dormancy", "Y" = "Dormancy")) +
  labs(y = "Count (Seed dormancy)", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(0, 50))

weightData <- conifer[, c("mastEvent", "logSeedWeight")]
weightData$mastEvent <- as.factor(weightData$mastEvent)

weight <- ggplot(weightData, aes(x = mastEvent, y = logSeedWeight, fill = mastEvent)) +
  geom_violin(trim = FALSE, alpha = 0.5, colour = "black") +
  geom_point(position = position_jitter(width = 0.08), size = 1.5, alpha = 0.7)  +
annotate("text", x = 1.5, y = 11, label = getAnnotation("Seed weight (log)", conifer_results, weightData),
         size = 2)+
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#A9D5B1", "1"="#ED562C"), 
    labels = c("0"="No", "1"="Yes")) + scale_x_discrete(labels = c("0" = "Non-masting", "1" = "Masting")) +
  labs(y = "Seed Weight (log)", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(NA, 13))

fruitData <- conifer[, c("mastEvent", "logFruit")]
fruitData$mastEvent <- as.factor(fruitData$mastEvent)

fruit <- ggplot(fruitData, aes(x = mastEvent, y = logFruit, fill = mastEvent)) +
  geom_violin(trim = FALSE, alpha = 0.5, colour = "black") + geom_point(position = position_jitter(width = 0.08), size = 1.5, alpha = 0.7) +
annotate("text", x = 1.5, y = 5, label = getAnnotation("Fruit size (log)", conifer_results, fruitData),
         size = 2)+
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#A9D5B1", "1"="#ED562C"), 
    labels = c("0"="No", "1"="Yes") 
  ) + scale_x_discrete(labels = c("0" = "Non-masting", "1" = "Masting")) +
  labs(y = "Fruit size (log)", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(NA, 6))

oilData <- conifer[, c("mastEvent", "oilContent")]
oilData$mastEvent <- as.factor(oilData$mastEvent)

oil <- ggplot(oilData, aes(x = mastEvent, y = oilContent, fill = mastEvent)) +
  geom_violin(trim = FALSE, alpha = 0.5, colour = "black") +
  geom_point(position = position_jitter(width = 0.08), size = 1.5, alpha = 0.7) +
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#A9D5B1", "1"="#ED562C"), 
    labels = c("0"="No", "1"="Yes") , guide = "none"
  ) + scale_x_discrete(labels = c("0" = "Non-masting", "1" = "Masting")) +
  labs(y = "Oil content %", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(NA, 100))

#pollData <- conifer[, c("mastEvent", "pollination")]
#pollData$mastEvent <- as.factor(pollData$mastEvent)
#pollData$pollination <- as.factor(pollData$pollination)
# Create ggplot object
#poll <- ggplot(pollData, aes(x = pollination, fill = mastEvent)) +
#  geom_bar(position = "dodge", colour = "black")  +
#  scale_fill_manual(
#    name = "Masting", 
#    values = c("0"="#A9D5B1", "1"="#ED562C"), 
#    labels = c("0"="No", "1"="Yes") 
#  ) + scale_x_discrete(labels = c("animals" = "Animal", "wind" = #"Wind", "wind and animals"="Both")) +
#  labs(y = "Count", x = NULL) +
 # theme_classic(base_size = 12) + scale_y_continuous(limits = c(0, 35#))

repData <- conifer[, c("mastEvent", "typeMonoOrDio")]
repData$mastEvent <- as.factor(repData$mastEvent)
repData$typeMonoOrDio <- as.factor(repData$typeMonoOrDio)
# Create ggplot object
rep <- ggplot(repData, aes(x = typeMonoOrDio, fill = mastEvent)) +
  geom_bar(position = "dodge", colour = "black") +
annotate("text", x = 1.5, y = 55, label = getAnnotation("Reproductive type", conifer_results, repData),
         size = 2)+
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#A9D5B1", "1"="#ED562C"), 
    labels = c("0"="No", "1"="Yes") 
  ) +
  labs(y = "Count (Reproductive type)", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(0, 60))

droughtData <- conifer[, c("mastEvent", "droughtTolerance")]
droughtData$mastEvent <- as.factor(droughtData$mastEvent)
droughtData$droughtTolerance <- as.factor(droughtData$droughtTolerance)
# Create ggplot object
drought <- ggplot(droughtData, aes(x = droughtTolerance, fill = mastEvent)) +
  geom_bar(position = "dodge", colour = "black")  +   annotate("text", x = 2.5, y = 23,label = getAnnotation("Drought tolerance", conifer_results, droughtData),  size = 2) +
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#A9D5B1", "1"="#ED562C"), 
    labels = c("0"="No", "1"="Yes") 
  ) +
  labs(y = "Count (Drought tolerance)", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(0, 25))

leafData <- conifer[, c("mastEvent", "leafLongevity")]
leafData$mastEvent <- as.factor(leafData$mastEvent)

leaf <- ggplot(leafData, aes(x = mastEvent, y = leafLongevity, fill = mastEvent)) +
  geom_violin(trim = FALSE, alpha = 0.5, colour = "black") +
  geom_point(position = position_jitter(width = 0.08), size = 1.5, alpha = 0.7)  +
  annotate("text", x = 1.5, y = 11, 
           label = getAnnotation("Leaf longevity (years)", conifer_results, leafData),
           size = 2) +
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#A9D5B1", "1"="#ED562C"), 
    labels = c("0"="No", "1"="Yes"))  + scale_x_discrete(labels = c("0" = "Non-masting", "1" = "Masting")) +
  labs(y = "Leaf longevity (year)", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(NA, 12))

predation <- weight + fruit + oil + disp + dorm 
pollination <- rep + plot_layout(ncol = 3, widths = c(1,1,1))+ guides(color = "none")+ guides(color = "none")
resource <- drought + leaf + plot_layout(ncol = 3, widths = c(1,1,1))

final <- predation / pollination / resource +
  plot_layout(heights = c(1, 0.5, 0.5), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(face = "bold", size = 8)
  )
ggsave(
  filename = "output/figures/rawDataWithStatsCon.pdf",
  plot = final,
  width = 12,
  height = 12
)
