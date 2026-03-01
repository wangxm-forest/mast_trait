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
# Extract key results from a glm object
tidy_glm <- function(model) {
  s <- summary(model)$coefficients
  out <- data.frame(
    term = rownames(s),
    estimate = s[,1],
    std_error = s[,2],
    z_value = s[,3],
    p_value = s[,4],
    N = nobs(model),
    row.names = NULL
  )
  return(out)
}

# Extract key results from a lm object
tidy_lm <- function(model) {
  s <- summary(model)$coefficients
  out <- data.frame(
    term = rownames(s),
    estimate = s[, 1],
    std_error = s[, 2],
    z_value = s[, 3],
    p_value = s[, 4],
    N = nobs(model),
    row.names = NULL
  )
  return(out)
}

# Run glm model and attach metadata
run_glm <- function(formula, data, method, trait_name) {
  model <- glm(formula, data = data, method = method)
  tbl <- tidy_glm(model)
  tbl$trait <- trait_name
  tbl$model_formula <- deparse(formula)
  return(tbl)
}

# Run phyloglm model and attach metadata
run_phyloglm <- function(formula, data, phy, method, trait_name) {
  model <- phyloglm(formula, phy = phy, data = data, method = method)
  tbl <- tidy_phyloglm(model)
  tbl$trait <- trait_name
  tbl$model_formula <- deparse(formula)
  return(tbl)
}

# Run LM model and attach metadata
run_lm <- function(formula, data, trait_name) {
  model <- lm(formula, data = data)
  tbl <- tidy_lm(model)
  tbl$trait <- trait_name
  tbl$model_formula <- deparse(formula)
  return(tbl)
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

table_grob <- tableGrob(
  conifer_results,
  rows = NULL,
  theme = ttheme_minimal(
    core = list(fg_params = list(fontsize = 9)),
    colhead = list(fg_params = list(fontsize = 10, fontface = "bold"))
  )
)
ggsave(
  filename = "output/phyloglmConifer.pdf",
  plot = table_grob,
  width = 8,
  height = 4
)

angio_results <- clean_results(angio_results)

table_grob <- tableGrob(
  angio_results,
  rows = NULL,
  theme = ttheme_minimal(
    core = list(fg_params = list(fontsize = 9)),
    colhead = list(fg_params = list(fontsize = 10, fontface = "bold"))
  )
)
ggsave(
  filename = "output/phyloglmAngio.pdf",
  plot = table_grob,
  width = 8,
  height = 4
)


### Use conifer and angiosperm as a fixed effect ----


# Model definitions
model_list <- list(
  list(name="Seed dispersal", formula=mastEvent ~ seedDispersal * group, data=d, method="brglmFit"),
 
  list(name="Pollination", formula=mastEvent ~ pollination * group, data=d, method="brglmFit"),
  
  list(name="Seed dormancy", formula=mastEvent ~ seedDormancy * group, data=d, method="brglmFit"),
 
  list(name="Mono/Dio", formula=mastEvent ~ typeMonoOrDio * group, data=d, method="brglmFit"),
 
  list(name="Seed weight (log)",    formula=mastEvent ~ logSeedWeight * group, data=d, method="brglmFit"),
  
  list(name="Fruit size (log)",     formula=mastEvent ~ logFruit * group, data=d, method="brglmFit"),
  
  list(name="Seed size (log)",      formula=mastEvent ~ logSeedSize * group, data=d, method="brglmFit"),
  
  list(name="Oil content %",      formula=mastEvent ~ oilContent * group, data=d, method="brglmFit"),
  
  list(name="Leaf longevity (years)", formula=mastEvent ~ leafLongevity * group, data=d, method="brglmFit"),
  
  list(name="Drought tolerance",    formula=mastEvent ~ droughtTolerance * group, data=d, method="brglmFit")
)

# Run models

results_glm <- NULL

for (m in model_list) {
  cat("Running model:", m$name, "\n")
  
  tbl <- run_glm(
    formula = m$formula,
    data    = m$data,
    method  = m$method,
    trait_name = m$name
  )
  
  results_glm <- rbind(results_glm, tbl)
}

write.csv(results_glm, "output/glmResults.csv", row.names = FALSE)

####Plot the probability for all traits####
#invlogit <- function(x) 1 / (1 + exp(-x))

my_colors <- c("angiosperm" = "#95B958",
               "conifer"    = "#6194BF")

custom_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.25),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.major = element_line(color = "gray85"),
    panel.grid.minor = element_blank()
  )

# Extract legend
get_legend <- function(myplot) {
  g <- ggplotGrob(myplot)
  leg <- g$grobs[which(sapply(g$grobs, function(x) x$name) == "guide-box")]
  leg[[1]]
}

plot_log_with_signif <- function(pred_df, trait_name, results_df, trait_var) {
  
  # Filter results for this trait (exclude intercept)
  trait_results <- results_df %>%
    filter(trait == trait_name & term != "(Intercept)")
  # Create significance column: default "ns"
  pred_df$signif <- "ns"
  
  # Assign "*" for significant contrasts (p < 0.05)
  for (i in 1:nrow(trait_results)) {
    term_name <- trait_results$term[i]
    p_val <- trait_results$p_value[i]
    
    # remove trait_var prefix to match trait_level in pred_df
    lvl <- gsub(trait_var, "", term_name)
    
    pred_df$signif[pred_df$trait_level == lvl] <- ifelse(p_val < 0.05, "*", "ns")
  }
  
  # Plot with significance as text above points
  ggplot(pred_df, aes(x = trait_level, y = logit, color = group, group = group)) +
    geom_point(size = 3) +
    geom_text(aes(label = signif), vjust = -0.5, size = 5, show.legend = FALSE,
              position = position_dodge(width = 1)) +
    labs(
      x = trait_name,
      y = "Log odds"
    ) + scale_color_manual(values = my_colors) +
    theme_classic()
}

# Pollination
poll <- subset(results_glm, trait == "Pollination")
coefs <- setNames(poll$estimate, poll$term)

predPoll <- expand.grid(
  trait_level = c("animal", "wind", "wind and animals"),
  group = c("angiosperm", "conifer")
)

predPoll$logit <- NA

predPoll$logit[predPoll$group == "angiosperm" &
                 predPoll$trait_level == "animal"] <- coefs["(Intercept)"]

predPoll$logit[predPoll$group == "angiosperm" &
                 predPoll$trait_level == "wind"] <-
  coefs["(Intercept)"] + coefs["pollinationwind"]

predPoll$logit[predPoll$group == "angiosperm" &
                 predPoll$trait_level == "wind and animals"] <-
  coefs["(Intercept)"] + coefs["pollinationwind and animals"]


# Seed dispersal
disp <- subset(results_glm, trait == "Seed dispersal")
coefs <- setNames(disp$estimate, disp$term)

predDisp <- expand.grid(
  trait_level = c("abiotic", "biotic", "both"),
  group = c("angiosperm", "conifer")
)

predDisp$logit <- NA

predDisp$logit[predDisp$group == "angiosperm" & predDisp$trait_level == "abiotic"] <- coefs["(Intercept)"]

# Angiosperm effects
predDisp$logit[predDisp$group == "angiosperm" & predDisp$trait_level == "biotic"] <- 
  coefs["(Intercept)"] + coefs["seedDispersalbiotic"]

predDisp$logit[predDisp$group == "angiosperm" & predDisp$trait_level == "both"] <- 
  coefs["(Intercept)"] + coefs["seedDispersalboth"]

# Conifer effects
predDisp$logit[predDisp$group == "conifer" & predDisp$trait_level == "abiotic"] <- 
  coefs["(Intercept)"] + coefs["groupconifer"]

predDisp$logit[predDisp$group == "conifer" & predDisp$trait_level == "biotic"] <- 
  coefs["(Intercept)"] + coefs["groupconifer"] +
  coefs["seedDispersalbiotic"] + coefs["seedDispersalbiotic:groupconifer"]

predDisp$logit[predDisp$group == "conifer" & predDisp$trait_level == "both"] <- 
  coefs["(Intercept)"] + coefs["groupconifer"] +
  coefs["seedDispersalboth"] + coefs["seedDispersalboth:groupconifer"]

# Seed dormancy
dorm <- subset(results_glm, trait == "Seed dormancy")
coefs <- setNames(dorm$estimate, dorm$term)

predDorm <- expand.grid(
  trait_level = c("N","Y"),
  group = c("angiosperm", "conifer")
)

predDorm$logit <- NA

predDorm$logit[predDorm$group == "angiosperm" &
                 predDorm$trait_level == "N"] <- coefs["(Intercept)"]
predDorm$logit[predDorm$group == "angiosperm" &
                 predDorm$trait_level == "Y"] <- coefs["(Intercept)"] + coefs["seedDormancyY"]
predDorm$logit[predDorm$group == "conifer" &
                 predDorm$trait_level == "N"] <- coefs["(Intercept)"] + coefs["groupconifer"]
predDorm$logit[predDorm$group == "conifer" &
                 predDorm$trait_level == "Y"] <- coefs["(Intercept)"] + coefs["seedDormancyY"] + coefs["groupconifer"] + coefs["seedDormancyY:groupconifer"]


# Mono/Dio
mono <- subset(results_glm, trait == "Mono/Dio")
coefs <- setNames(mono$estimate, mono$term)

predMono <- expand.grid(
  trait_level = c("Dioecious", "Monoecious", "Polygamous"),
  group = c("angiosperm", "conifer")
)

predMono$logit <- NA
predMono$logit[predMono$group == "angiosperm" &
                 predMono$trait_level == "Dioecious"] <- coefs["(Intercept)"]
predMono$logit[predMono$group == "angiosperm" &
                 predMono$trait_level == "Monoecious"] <- coefs["(Intercept)"] + coefs["typeMonoOrDioMonoecious"]
predMono$logit[predMono$group == "angiosperm" &
                 predMono$trait_level == "Polygamous"] <- coefs["(Intercept)"] + coefs["typeMonoOrDioPolygamous"]
predMono$logit[predMono$group == "conifer" &
                 predMono$trait_level == "Dioecious"] <- coefs["(Intercept)"] + coefs["groupconifer"]
predMono$logit[predMono$group == "conifer" &
                 predMono$trait_level == "Monoecious"] <- coefs["(Intercept)"] + coefs["groupconifer"] + coefs["typeMonoOrDioMonoecious:groupconifer"]

# Drought tolerance
drought <- subset(results_glm, trait == "Drought tolerance")
coefs <- setNames(drought$estimate, drought$term)

preddrought <- expand.grid(
  trait_level = c("High", "Low", "Moderate"),
  group = c("angiosperm", "conifer")
)

preddrought$logit <- NA
preddrought$logit[preddrought$group == "angiosperm" &
                    preddrought$trait_level == "High"] <- coefs["(Intercept)"]
preddrought$logit[preddrought$group == "angiosperm" &
                    preddrought$trait_level == "Low"] <- coefs["(Intercept)"] + coefs["droughtToleranceLow"]
preddrought$logit[preddrought$group == "angiosperm" &
                    preddrought$trait_level == "Moderate"] <- coefs["(Intercept)"] + coefs["droughtToleranceModerate"]
preddrought$logit[preddrought$group == "conifer" &
                    preddrought$trait_level == "High"] <- coefs["(Intercept)"] + coefs["groupconifer"]
preddrought$logit[preddrought$group == "conifer" &
                    preddrought$trait_level == "Low"] <- coefs["(Intercept)"] + coefs["groupconifer"] + coefs["droughtToleranceLow:groupconifer"]
preddrought$logit[preddrought$group == "conifer" &
                    preddrought$trait_level == "Moderate"] <- coefs["(Intercept)"] + coefs["groupconifer"] + coefs["droughtToleranceModerate:groupconifer"]

#Plotting
plotPoll <- plot_log_with_signif(predPoll, 
                      trait_name = "Pollination", 
                      results_df = results_glm, 
                      trait_var = "pollination") + theme(legend.position = "none")

plotDisp <- plot_log_with_signif(predDisp, 
                      trait_name = "Seed dispersal", 
                      results_df = results_glm, 
                      trait_var = "seedDispersal") + theme(legend.position = "none",axis.title.y = element_blank()) 

plotDorm <- plot_log_with_signif(predDorm, 
                      trait_name = "Seed dormancy", 
                      results_df = results_glm, 
                      trait_var = "seedDormancy") + theme(legend.position = "none",axis.title.y = element_blank())

plotMono <- plot_log_with_signif(predMono, 
                      trait_name = "Mono/Dio", 
                      results_df = results_glm, 
                      trait_var = "typeMonoOrDio") + theme(legend.position = "none",axis.title.y = element_blank())

plotDrought <- plot_log_with_signif(preddrought, 
                      trait_name = "Drought tolerance", 
                      results_df = results_glm, 
                      trait_var = "droughtTolerance")
shared_legend <- get_legend(plotDrought)
plotDrought <- plotDrought + theme(legend.position = "none",axis.title.y = element_blank())

pdf("output/figures/glmCat.pdf", width = 20, height = 5)
grid.arrange(plotPoll, plotDisp, plotDorm, plotMono,plotDrought, shared_legend, ncol = 6)
dev.off()
## Continuous traits
#Seed weight
mean <- mean(d$logSeedWeight, na.rm =TRUE)
seedweight <- subset(results_glm, trait == "Seed weight (log)")
coefs <- setNames(seedweight$estimate, seedweight$term)

seedweight <- expand.grid(
  mean = mean,
  group = c("angiosperm", "conifer")
)

seedweight$logit <- with(seedweight,
                   coefs["(Intercept)"] +
                     coefs["logSeedWeight"] * mean +
                     ifelse(group == "conifer", coefs["groupconifer"], 0) +
                     ifelse(group == "conifer", coefs["logSeedWeight:groupconifer"] * mean, 0)
)

#Seed size
mean <- mean(d$logSeedSize, na.rm =TRUE)
seedsize <- subset(results_glm, trait == "Seed size (log)")
coefs <- setNames(seedsize$estimate, seedsize$term)

seedsize <- expand.grid(
  mean = mean,
  group = c("angiosperm", "conifer")
)

seedsize$logit <- with(seedsize,
                   coefs["(Intercept)"] +
                     coefs["logSeedSize"] * mean +
                     ifelse(group == "conifer", coefs["groupconifer"], 0) +
                     ifelse(group == "conifer", coefs["logSeedSize:groupconifer"] * mean, 0)
)


#Fruit size
mean <- mean(d$logFruit, na.rm =TRUE)
fruitsize <- subset(results_glm, trait == "Fruit size (log)")
coefs <- setNames(fruitsize$estimate, fruitsize$term)

fruitsize <- expand.grid(
  mean = mean,
  group = c("angiosperm", "conifer")
)

fruitsize$logit <- with(fruitsize,
                       coefs["(Intercept)"] +
                         coefs["logFruit"] * mean +
                         ifelse(group == "conifer", coefs["groupconifer"], 0) +
                         ifelse(group == "conifer", coefs["logFruit:groupconifer"] * mean, 0)
)


#Oil content
mean <- mean(d$oilContent, na.rm =TRUE)
oilcontent <- subset(results_glm, trait == "Oil content %")
coefs <- setNames(oilcontent$estimate, oilcontent$term)

oilcontent <- expand.grid(
  mean = mean,
  group = c("angiosperm", "conifer")
)

oilcontent$logit <- with(oilcontent,
                        coefs["(Intercept)"] +
                          coefs["oilContent"] * mean +
                          ifelse(group == "conifer", coefs["groupconifer"], 0) +
                          ifelse(group == "conifer", coefs["oilContent:groupconifer"] * mean, 0)
)

#leaf longevity
mean <- mean(d$leafLongevity, na.rm =TRUE)
leaflongevity <- subset(results_glm, trait == "Leaf longevity (years)")
coefs <- setNames(leaflongevity$estimate, leaflongevity$term)

leaflongevity <- expand.grid(
  mean = mean,
  group = c("angiosperm", "conifer")
)

leaflongevity$logit <- with(leaflongevity,
                         coefs["(Intercept)"] +
                           coefs["leafLongevity"] * mean +
                           ifelse(group == "conifer", coefs["groupconifer"], 0) +
                           ifelse(group == "conifer", coefs["leafLongevity:groupconifer"] * mean, 0)
)

plot_continuous_traits_combined <- function(pred_df) {
  
  ggplot(pred_df, aes(x = trait_name, y = logit, color = group, group = group)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_text(aes(label = signif),
              position = position_dodge(width = 0.5),
              vjust = -0.8,
              size = 5,
              show.legend = FALSE) +
    labs(
      x = "Trait",
      y = "Log odds"
    ) +
    scale_color_manual(values = my_colors) +
    theme_classic()
}

# Add trait_name column to each pred_df
seedweight$trait_name <- "Seed weight (log)"
seedsize$trait_name <- "Seed size (log)"
fruitsize$trait_name <- "Fruit size (log)"
oilcontent$trait_name <- "Oil content %"
leaflongevity$trait_name <- "Leaf longevity (years)"

# Combine into one dataframe
pred_all <- bind_rows(seedweight, seedsize, fruitsize, oilcontent, leaflongevity)
pred_all$signif <- c("*","*","ns","ns","ns","ns","ns","ns","ns","ns")

plot_combined <- plot_continuous_traits_combined(pred_all)
plot_combined

pdf("output/figures/glmCon.pdf", width = 10, height = 10)
grid.arrange(plotSeedweight, plotSeedsize, plotFruitsize, plotOilcontent, plotLeaflongevity, shared_legend, ncol = 6)
dev.off()

### Plot the mean and SE ----
plot_df <- d

all_df <- d
all_df$group <- "All species"

plot_df <- rbind(plot_df, all_df)

# Fruit size
mean_vals <- tapply(plot_df$logFruitStd, plot_df$group, mean, na.rm = TRUE)
sd_vals   <- tapply(plot_df$logFruitStd, plot_df$group, sd,   na.rm = TRUE)
n_vals    <- tapply(plot_df$logFruitStd, plot_df$group, function(x) sum(!is.na(x)))

se_vals <- sd_vals / sqrt(n_vals)

summary_df <- data.frame(
  group = names(mean_vals),
  mean  = as.numeric(mean_vals),
  se    = as.numeric(se_vals),
  n     = as.numeric(n_vals)
)

# control order on x-axis
summary_df$group <- factor(
  summary_df$group,
  levels = c("angiosperm", "conifer", "All species")
)

fruit_size <- ggplot(summary_df, aes(x = group, y = mean, color = group)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.05
  ) +
  geom_text(
    aes(label = paste0("n = ", n)),
    nudge_x = 0.2,
    size = 5,
    show.legend = FALSE
  ) +
  custom_theme + scale_color_manual(values = my_colors) +
  labs(
    x = "",
    y = "Fruit size (log)"
  ) + theme(legend.position = "none")

fruit_size_raw <- ggplot(plot_df, aes(x = group, y = logFruitStd, color = group)) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15) +
  scale_color_manual(values = my_colors) +
  theme_classic() +
  labs(
    x = "",
    y = "Fruit size (log)"
  ) + theme(legend.position = "none")
# Seed size
mean_vals <- tapply(plot_df$logSeedSizeStd, plot_df$group, mean, na.rm = TRUE)
sd_vals   <- tapply(plot_df$logSeedSizeStd, plot_df$group, sd,   na.rm = TRUE)
n_vals    <- tapply(plot_df$logSeedSizeStd, plot_df$group, function(x) sum(!is.na(x)))

se_vals <- sd_vals / sqrt(n_vals)

summary_df <- data.frame(
  group = names(mean_vals),
  mean  = as.numeric(mean_vals),
  se    = as.numeric(se_vals),
  n     = as.numeric(n_vals)
)

# control order on x-axis
summary_df$group <- factor(
  summary_df$group,
  levels = c("angiosperm", "conifer", "All species")
)

seed_size <- ggplot(summary_df, aes(x = group, y = mean, color = group)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.05
  ) +
  geom_text(
    aes(label = paste0("n = ", n)),
    nudge_x = 0.2,
    size = 5,
    show.legend = FALSE
  ) +
  custom_theme + scale_color_manual(values = my_colors) +
  labs(
    x = "",
    y = "Seed size (log)"
  ) + theme(legend.position = "none")
seed_size_raw <- ggplot(plot_df, aes(x = group, y = logSeedSizeStd, color = group)) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15) +
  scale_color_manual(values = my_colors) +
  theme_classic() +
  labs(
    x = "",
    y = "Seed size (log)"
  ) + theme(legend.position = "none")
# Seed weight
mean_vals <- tapply(plot_df$logSeedWeightStd, plot_df$group, mean, na.rm = TRUE)
sd_vals   <- tapply(plot_df$logSeedWeightStd, plot_df$group, sd,   na.rm = TRUE)
n_vals    <- tapply(plot_df$logSeedWeightStd, plot_df$group, function(x) sum(!is.na(x)))

se_vals <- sd_vals / sqrt(n_vals)

summary_df <- data.frame(
  group = names(mean_vals),
  mean  = as.numeric(mean_vals),
  se    = as.numeric(se_vals),
  n     = as.numeric(n_vals)
)

# control order on x-axis
summary_df$group <- factor(
  summary_df$group,
  levels = c("angiosperm", "conifer", "All species")
)

seed_weight <- ggplot(summary_df, aes(x = group, y = mean, color = group)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.05
  ) +
  geom_text(
    aes(label = paste0("n = ", n)),
    nudge_x = 0.2,
    size = 5,
    show.legend = FALSE
  ) +
  custom_theme + scale_color_manual(values = my_colors) +
  labs(
    x = "",
    y = "Seed weight (log)"
  ) + theme(legend.position = "none")
seed_weight_raw <- ggplot(plot_df, aes(x = group, y = logSeedWeightStd, color = group)) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15) +
  scale_color_manual(values = my_colors) +
  theme_classic() +
  labs(
    x = "",
    y = "Seed weight (log)"
  ) + theme(legend.position = "none")

# oil content
plot_df <- d

all_df <- d
all_df$group <- "All species"

plot_df <- rbind(plot_df, all_df)
mean_vals <- tapply(plot_df$oilContent, plot_df$group, mean, na.rm = TRUE)
sd_vals   <- tapply(plot_df$oilContent, plot_df$group, sd,   na.rm = TRUE)
n_vals    <- tapply(plot_df$oilContent, plot_df$group, function(x) sum(!is.na(x)))

se_vals <- sd_vals / sqrt(n_vals)

summary_df <- data.frame(
  group = names(mean_vals),
  mean  = as.numeric(mean_vals),
  se    = as.numeric(se_vals),
  n     = as.numeric(n_vals)
) 

# control order on x-axis
summary_df$group <- factor(
  summary_df$group,
  levels = c("angiosperm", "conifer", "All species")
)

oil_content <- ggplot(summary_df, aes(x = group, y = mean, color = group)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.05
  ) +
  geom_text(
    aes(label = paste0("n = ", n)),
    nudge_x = 0.2,
    size = 5,
    show.legend = FALSE
  ) +
  custom_theme + scale_color_manual(values = my_colors) + 
  labs(
    x = "",
    y = "Oil content %"
  ) + theme(legend.position = "none")
oil_content_raw <- ggplot(plot_df, aes(x = group, y = oilContent, color = group)) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15) +
  scale_color_manual(values = my_colors) +
  theme_classic() +
  labs(
    x = "",
    y = "Oil content %"
  ) + theme(legend.position = "none")
# leaf longevity
plot_df <- d
mean_vals <- tapply(plot_df$leafLongevity, plot_df$group, mean, na.rm = TRUE)
sd_vals   <- tapply(plot_df$leafLongevity, plot_df$group, sd,   na.rm = TRUE)
n_vals    <- tapply(plot_df$leafLongevity, plot_df$group, function(x) sum(!is.na(x)))

se_vals <- sd_vals / sqrt(n_vals)

summary_df <- data.frame(
  group = names(mean_vals),
  mean  = as.numeric(mean_vals),
  se    = as.numeric(se_vals),
  n     = as.numeric(n_vals)
)

# control order on x-axis
summary_df$group <- factor(
  summary_df$group,
  levels = c("angiosperm", "conifer")
)

leaf_longevity <- ggplot(summary_df, aes(x = group, y = mean, color = group)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.05
  ) +
  geom_text(
    aes(label = paste0("n = ", n)),
    nudge_x = 0.2,
    size = 5,
    show.legend = FALSE
  ) +
  custom_theme + scale_color_manual(values = my_colors) + 
  labs(
    x = "",
    y = "Leaf longevity (year)"
  ) + theme(legend.position = "none")
leaf_longevity_raw <- ggplot(plot_df, aes(x = group, y = leafLongevity, color = group)) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15) +
  scale_color_manual(values = my_colors) +
  theme_classic() +
  labs(
    x = "",
    y = "Leaf longevity (year)"
  ) + theme(legend.position = "none")
pdf("output/figures/meanSE.pdf", width = 25, height = 10)
grid.arrange(fruit_size, seed_size, seed_weight, oil_content,leaf_longevity, nrow = 2, ncol = 3)
dev.off()
pdf("output/figures/meanSERaw.pdf", width = 25, height = 10)
grid.arrange(fruit_size_raw, seed_size_raw, seed_weight_raw, oil_content_raw,leaf_longevity_raw, nrow = 2, ncol = 3)
dev.off()

#pgls
## Make a function to subset traits
subsetTrait <- function(data, cols) {
  dataSub <- data[, cols, drop = FALSE]
  dataSub <- dataSub[complete.cases(dataSub), ]
  dataSub
}

## Make a function to run the model on different traits
runPgls <- function(data, phy, traits,
                         response = "mastCycleAveLog",
                         speciesCol = "latbi",
                         minN = 5) {
  
  models <- list()
  
  for (tr in traits) {message("Running PGLS for trait: ", tr)
    
    cols <- c(speciesCol, response, tr)
    dSub <- subsetTrait(data, cols)
    
    if (nrow(dSub) < minN) {
      message("  Skipped (too few species)")
      next
    }
    
    ## Comparative data
    comp <- comparative.data(
      phy = phy,
      data = dSub,
      names.col = "latbi",
      vcv = TRUE,
      warn.dropped = TRUE
    )
    
    ## Fit PGLS
    form <- as.formula(paste(response, "~", tr))
    models[[tr]] <- pgls(form, data = comp, lambda = "ML")
  }
  
  models
}

## Extract results
extractResults <- function(models) {
  
  do.call(
    rbind,
    lapply(names(models), function(tr) {
      
      m <- models[[tr]]
      coef <- as.data.frame(summary(m)$coefficients)
      
      data.frame(
        trait     = tr,
        term      = rownames(coef),
        estimate  = coef$Estimate,
        std_error = coef$`Std. Error`,
        p_value   = coef$`Pr(>|t|)`,
        lambda    = m$param["lambda"],
        N         = m$n,
        row.names = NULL
      )
    })
    
    
  )
}

# Data prep
##response has to be Gaussian distribution
conifer$mastCycleAveLog <- log(conifer$mastCycleAve)
angio$mastCycleAveLog <- log(angio$mastCycleAve)

##traits of interest
conTraits <- c(
  "seedDispersal",
  "seedDormancy",
  "typeMonoOrDio",
  "droughtTolerance",
  "leafLongevity",
  "oilContent",
  "logFruitStd",
  "logSeedSizeStd",
  "logSeedWeightStd"
)

angioTraits <- c(
  "seedDispersal",
  "seedDormancy",
  "typeMonoOrDio",
  "droughtTolerance",
  "pollination",
  "leafLongevity",
  "oilContent",
  "logFruitStd",
  "logSeedSizeStd",
  "logSeedWeightStd"
)

# Run pgls
pglsConifer <- runPgls(
  data = conifer,
  phy = phyconifer,
  traits = conTraits
)

pglsAngio <- runPgls(
  data = angio,
  phy = phyangio,
  traits = angioTraits
)

# Extract results
resultsConifer <- extractResults(pglsConifer)
resultsAngio <- extractResults(pglsAngio)

# Clean the results
numCols <- c("estimate", "std_error", "p_value", "lambda")

resultsConifer[numCols] <- lapply(resultsConifer[numCols], round, 3)
resultsAngio[numCols] <- lapply(resultsAngio[numCols], round, 3)

traitLabels <- c(
  seedDispersal = "Dispersal mode",
  seedDormancy = "Seed dormancy",
  typeMonoOrDio = "Reproductive type",
  droughtTolerance = "Drought tolerance",
  pollination = "Pollination mode",
  leafLongevity = "Leaf longevity (years)",
  oilContent = "Oil content (%)",
  logSeedWeightStd = "Seed weight (log)",
  logFruitStd = "Fruit size (log)",
  logSeedSizeStd = "Seed size (log)"
)

termLabels <- c(
  "(Intercept)" = "Intercept",
  "seedDispersalbiotic" = "Biotic vs abiotic",
  "seedDispersalboth" = "Both vs abiotic",
  "seedDormancyY" = "Dormant vs non-dormant",
  "typeMonoOrDioMonoecious" = "Monoecious vs dioecious",
  "typeMonoOrDioPolygamous" = "Polygamous vs dioecious",
  "droughtToleranceLow" = "Low vs high drought tolerance",
  "droughtToleranceModerate" = "Moderate vs high drought tolerance",
  "pollinationwind" = "Wind vs animal pollination",
  "pollinationwind and animals" = "Wind+animal vs animal pollination",
  "leafLongevity" = "Leaf longevity (years)",
  "oilContent" = "Oil content (%)",
  "logFruitStd" = "Fruit size (log)",
  "logSeedSizeStd" = "Seed size (log)",
  "logSeedWeightStd" = "Seed weight (log)"
)

resultsConifer$trait <- traitLabels[resultsConifer$trait]
resultsConifer$term <- termLabels[resultsConifer$term]

resultsAngio$trait <- traitLabels[resultsAngio$trait]
resultsAngio$term <- termLabels[resultsAngio$term]

resultsConifer <- subset(resultsConifer, term != "Intercept")
resultsAngio <- subset(resultsAngio, term != "Intercept")

colnames(resultsConifer) <- c(
  "Trait", "Predictor", "Estimate", "Std_Error", "P_value", "Lambda", "N"
)
colnames(resultsAngio) <- c(
  "Trait", "Predictor", "Estimate", "Std_Error", "P_value", "Lambda", "N"
)

table_grob <- tableGrob(
  resultsConifer,
  rows = NULL,
  theme = ttheme_minimal(
    core = list(fg_params = list(fontsize = 9)),
    colhead = list(fg_params = list(fontsize = 10, fontface = "bold"))
  )
)
ggsave(
  filename = "output/pglsConifer.pdf",
  plot = table_grob,
  width = 8,
  height = 4
)

table_grob <- tableGrob(
  resultsAngio,
  rows = NULL,
  theme = ttheme_minimal(
    core = list(fg_params = list(fontsize = 9)),
    colhead = list(fg_params = list(fontsize = 10, fontface = "bold"))
  )
)
ggsave(
  filename = "output/pglsAngio.pdf",
  plot = table_grob,
  width = 8,
  height = 4
)

### including conifer/angiosperm as a fixed effect in a lm ----

d$mastCycleAveLog <- log(d$mastCycleAve)
# Model definitions
model_list <- list(
  list(name="Dispersal mode", formula=mastCycleAveLog ~ seedDispersal + group, data=d),
  
  list(name="Pollination mode", formula=mastCycleAveLog ~ pollination + group, data=d),
  
  list(name="Seed dormancy", formula=mastCycleAveLog ~ seedDormancy + group, data=d),
  
  list(name="Reproductive type", formula=mastCycleAveLog ~ typeMonoOrDio + group, data=d),
  
  list(name="Seed weight (log)",    formula=mastCycleAveLog ~ logSeedWeightStd + group, data=d),
  
  list(name="Fruit size (log)",     formula=mastCycleAveLog ~ logFruitStd + group, data=d),
  
  list(name="Seed size (log)",      formula=mastCycleAveLog ~ logSeedSizeStd + group, data=d),
  
  list(name="Oil content (%)",      formula=mastCycleAveLog ~ oilContent + group, data=d),
  
  list(name="Leaf longevity (years)", formula=mastCycleAveLog ~ leafLongevity + group, data=d),
  
  list(name="Drought tolerance",    formula=mastCycleAveLog ~ droughtTolerance + group, data=d)
)

# Run models

results_lm <- NULL

for (m in model_list) {
  cat("Running model:", m$name, "\n")
  
  tbl <- run_lm(
    formula = m$formula,
    data    = m$data,
    trait_name = m$name
  )
  
  results_lm <- rbind(results_lm, tbl)
}

#Make a table to present the results:
results_no_int <- results_lm[results_lm$term != "(Intercept)", ]

##results_no_int$signif <- cut(results_no_int$p_value, breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf), labels = c("***", "**", "*", ".", ""))

results_no_int$estimate  <- round(results_no_int$estimate,  3)
results_no_int$std_error <- round(results_no_int$std_error, 3)
results_no_int$t_value   <- round(results_no_int$t_value,   2)
results_no_int$p_value   <- round(results_no_int$p_value,  3)

final_table <- results_no_int[, c(
  "trait",
  "term",
  "estimate",
  "std_error",
  "p_value",
  "N"
)]

colnames(final_table) <- c(
  "Trait",
  "Predictor",
  "Estimate",
  "SE",
  "P",
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
final_table$Predictor <- gsub("groupconifer",  "Conifer compared to angiosperm", final_table$Predictor)

table_grob <- tableGrob(
  final_table,
  rows = NULL,
  theme = ttheme_minimal(
    core = list(fg_params = list(fontsize = 9)),
    colhead = list(fg_params = list(fontsize = 10, fontface = "bold"))
  )
)
ggsave(
  filename = "output/lmResults.pdf",
  plot = table_grob,
  width = 8,
  height = 8
)

###Plot the raw data with the phyloglm results ----
## Angiosperm
getAnnotation <- function(trait_name, model_results) {
  
  df <- model_results[model_results$Trait == trait_name, ]
  
  # Always keep alpha
  alpha_text <- paste0("Phylo α = ", unique(df$`Phylo α`))
  
  # Keep only predictors with P < 0.5
  df_sig <- df[df$P < 0.05, ]
  
  # If none meet threshold → return only alpha
  if(nrow(df_sig) == 0){
    return(alpha_text)
  }
  
  # If multiple predictors (categorical trait)
  if(length(unique(df$Predictor)) > 1){
    
    pred_text <- paste(
      df_sig$Predictor,
      ": Estimate=", round(df_sig$Estimate, 2),
      " ± ", round(df_sig$SE, 2),
      " P=", round(df_sig$P, 2),
      sep = "",
      collapse = "\n"
    )
    
    label <- paste(alpha_text, pred_text, sep = "\n")
    
  } else {
    
    # Single predictor
    pred_text <- paste0(
      df_sig$Predictor,
      ": Estimate=", round(df_sig$Estimate, 2),
      " ± ", round(df_sig$SE, 2),
      " P=", round(df_sig$P, 2)
    )
    
    label <- paste(alpha_text, pred_text, sep = "\n")
  }
  
  return(label)
}

getAnnotation("Dispersal mode", angio_results)
# "n=50, P=0.03,0.12"

dispData <- angio[, c("mastEvent", "seedDispersal")]
dispData$mastEvent <- as.factor(dispData$mastEvent)
dispData$seedDispersal <- as.factor(dispData$seedDispersal)
# Create ggplot object
disp <- ggplot(dispData, aes(x = seedDispersal, fill = mastEvent)) +
  geom_bar(position = "dodge", colour = "black") +
  annotate("text", x = 2.5, y = c(27,28), 
           label = c( "Nmast = 51, Nnon = 47", getAnnotation("Dispersal mode", angio_results)), size = 2) +
  labs(y = "Count", x = NULL) +
  theme_classic(base_size = 12) + scale_x_discrete(labels = c("abiotic" = "Abiotic", "biotic" = "Biotic", "both"="Both")) + scale_fill_manual(
    name = "Masting", 
    values = c("0"="#387E46", "1"="#C43142"), 
    labels = c("0"="No", "1"="Yes") 
  ) + scale_y_continuous(limits = c(0, 30))

dormData <- angio[, c("mastEvent", "seedDormancy")]
dormData$mastEvent <- as.factor(dormData$mastEvent)
dormData$seedDormancy <- as.factor(dormData$seedDormancy)
# Create ggplot object
dorm <- ggplot(dormData, aes(x = seedDormancy, fill = mastEvent)) +
  geom_bar(position = "dodge", colour = "black") +
  annotate("text", x = 2, y = 42, 
           label = getAnnotation("Seed dormancy", angio_results),
           size = 2) +
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#387E46", "1"="#C43142"), 
    labels = c("0"="No", "1"="Yes") 
  ) + scale_x_discrete(labels = c("N" = "No dormancy", "Y" = "Dormancy")) +
  labs(y = "Count", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(0, 45))

weightData <- angio[, c("mastEvent", "logSeedWeight")]
weightData$mastEvent <- as.factor(weightData$mastEvent)

weight <- ggplot(weightData, aes(x = mastEvent, y = logSeedWeight, fill = mastEvent)) +
  geom_violin(trim = FALSE, alpha = 0.5, colour = "black") +
  geom_point(position = position_jitter(width = 0.08), size = 1.5, alpha = 0.7)  +
  annotate("text", x = 1.5, y = 16, 
           label = getAnnotation("Seed weight (log)", angio_results),
           size = 2) +
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#387E46", "1"="#C43142"), 
    labels = c("0"="No", "1"="Yes")) + scale_x_discrete(labels = c("0" = "Non-masting", "1" = "Masting")) +
  labs(y = "Seed Weight (log)", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(NA, 18))

fruitData <- angio[, c("mastEvent", "logFruit")]
fruitData$mastEvent <- as.factor(fruitData$mastEvent)

fruit <- ggplot(fruitData, aes(x = mastEvent, y = logFruit, fill = mastEvent)) +
  geom_violin(trim = FALSE, alpha = 0.5, colour = "black") +
  annotate("text", x = 1.5, y = 5.5, 
           label = getAnnotation("Fruit size (log)", angio_results),
           size = 2) +
  geom_point(position = position_jitter(width = 0.08), size = 1.5, alpha = 0.7) +
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#387E46", "1"="#C43142"), 
    labels = c("0"="No", "1"="Yes") 
  ) + scale_x_discrete(labels = c("0" = "Non-masting", "1" = "Masting")) +
  labs(y = "Fruit size (log)", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(NA, 6))

oilData <- angio[, c("mastEvent", "oilContent")]
oilData$mastEvent <- as.factor(oilData$mastEvent)

oil <- ggplot(oilData, aes(x = mastEvent, y = oilContent, fill = mastEvent)) +
  geom_violin(trim = FALSE, alpha = 0.5, colour = "black") +
  annotate("text", x = 1.5, y = 110, 
           label = getAnnotation("Oil content (%)", angio_results),
           size = 2) +
  geom_point(position = position_jitter(width = 0.08), size = 1.5, alpha = 0.7) +
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#387E46", "1"="#C43142"), 
    labels = c("0"="No", "1"="Yes") 
  ) + scale_x_discrete(labels = c("0" = "Non-masting", "1" = "Masting")) +
  labs(y = "Oil content %", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(NA, 120))

pollData <- angio[, c("mastEvent", "pollination")]
pollData$mastEvent <- as.factor(pollData$mastEvent)
pollData$pollination <- as.factor(pollData$pollination)
# Create ggplot object
poll <- ggplot(pollData, aes(x = pollination, fill = mastEvent)) +
  geom_bar(position = "dodge", colour = "black")  +
  annotate("text", x = 2.5, y = 31, 
           label = getAnnotation("Pollination mode", angio_results),
           size = 2) +
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#387E46", "1"="#C43142"), 
    labels = c("0"="No", "1"="Yes") 
  ) + scale_x_discrete(labels = c("animals" = "Animal", "wind" = "Wind", "wind and animals"="Both")) +
  labs(y = "Count", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(0, 35))

repData <- angio[, c("mastEvent", "typeMonoOrDio")]
repData$mastEvent <- as.factor(repData$mastEvent)
repData$typeMonoOrDio <- as.factor(repData$typeMonoOrDio)
# Create ggplot object
rep <- ggplot(repData, aes(x = typeMonoOrDio, fill = mastEvent)) +
  geom_bar(position = "dodge", colour = "black") +
  annotate("text", x = 2.5, y = 46, 
           label = getAnnotation("Reproductive type", angio_results),
           size = 2) +
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#387E46", "1"="#C43142"), 
    labels = c("0"="No", "1"="Yes") 
  ) +
  labs(y = "Count", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(0, 50))

droughtData <- angio[, c("mastEvent", "droughtTolerance")]
droughtData$mastEvent <- as.factor(droughtData$mastEvent)
droughtData$droughtTolerance <- as.factor(droughtData$droughtTolerance)
# Create ggplot object
drought <- ggplot(droughtData, aes(x = droughtTolerance, fill = mastEvent)) +
  geom_bar(position = "dodge", colour = "black") +
  annotate("text", x = 2.5, y = 23, 
           label = getAnnotation("Drought tolerance", angio_results),
           size = 2) +
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#387E46", "1"="#C43142"), 
    labels = c("0"="No", "1"="Yes") 
  ) +
  labs(y = "Count", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(0, 25))

leafData <- angio[, c("mastEvent", "leafLongevity")]
leafData$mastEvent <- as.factor(leafData$mastEvent)

leaf <- ggplot(leafData, aes(x = mastEvent, y = leafLongevity, fill = mastEvent)) +
  geom_violin(trim = FALSE, alpha = 0.5, colour = "black") +
  geom_point(position = position_jitter(width = 0.08), size = 1.5, alpha = 0.7)  +
  annotate("text", x = 1.5, y = 3.5, 
           label = getAnnotation("Leaf longevity (years)", angio_results),
           size = 2) +
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#387E46", "1"="#C43142"), 
    labels = c("0"="No", "1"="Yes"))  + scale_x_discrete(labels = c("0" = "Non-masting", "1" = "Masting")) +
  labs(y = "Leaf longevity (year)", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(NA, 4))

predation <- weight + fruit + oil + disp + dorm
pollination <- poll + rep + plot_layout(ncol = 3, widths = c(1,1,1))
resource <- drought + leaf + plot_layout(ncol = 3, widths = c(1,1,1))

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
  annotate("text", x = 2.5, y = 23, 
           label = getAnnotation("Seed dispersal", conifer_results),
           size = 2) +
  labs(y = "Count", x = NULL) +
  theme_classic(base_size = 12) + scale_x_discrete(labels = c("abiotic" = "Abiotic", "biotic" = "Biotic", "both"="Both")) + scale_fill_manual(
    name = "Masting", 
    values = c("0"="#387E46", "1"="#C43142"), 
    labels = c("0"="No", "1"="Yes") 
  ) + scale_y_continuous(limits = c(0, 25))

dormData <- conifer[, c("mastEvent", "seedDormancy")]
dormData$mastEvent <- as.factor(dormData$mastEvent)
dormData$seedDormancy <- as.factor(dormData$seedDormancy)
# Create ggplot object
dorm <- ggplot(dormData, aes(x = seedDormancy, fill = mastEvent)) +
  geom_bar(position = "dodge", colour = "black") +
  annotate("text", x = 2, y = 45, 
           label = getAnnotation("Seed dormancy", conifer_results),
           size = 2) +
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#387E46", "1"="#C43142"), 
    labels = c("0"="No", "1"="Yes") 
  ) + scale_x_discrete(labels = c("N" = "No dormancy", "Y" = "Dormancy")) +
  labs(y = "Count", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(0, 50))

weightData <- conifer[, c("mastEvent", "logSeedWeight")]
weightData$mastEvent <- as.factor(weightData$mastEvent)

weight <- ggplot(weightData, aes(x = mastEvent, y = logSeedWeight, fill = mastEvent)) +
  geom_violin(trim = FALSE, alpha = 0.5, colour = "black") +
  geom_point(position = position_jitter(width = 0.08), size = 1.5, alpha = 0.7)  +
  annotate("text", x = 1.5, y = 11, 
           label = getAnnotation("Seed weight (log)", conifer_results),
           size = 2) +
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#387E46", "1"="#C43142"), 
    labels = c("0"="No", "1"="Yes")) + scale_x_discrete(labels = c("0" = "Non-masting", "1" = "Masting")) +
  labs(y = "Seed Weight (log)", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(NA, 13))

fruitData <- conifer[, c("mastEvent", "logFruit")]
fruitData$mastEvent <- as.factor(fruitData$mastEvent)

fruit <- ggplot(fruitData, aes(x = mastEvent, y = logFruit, fill = mastEvent)) +
  geom_violin(trim = FALSE, alpha = 0.5, colour = "black") +
  annotate("text", x = 1.5, y = 5, 
           label = getAnnotation("Fruit size (log)", conifer_results),
           size = 2) +
  geom_point(position = position_jitter(width = 0.08), size = 1.5, alpha = 0.7) +
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#387E46", "1"="#C43142"), 
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
    values = c("0"="#387E46", "1"="#C43142"), 
    labels = c("0"="No", "1"="Yes") 
  ) + scale_x_discrete(labels = c("0" = "Non-masting", "1" = "Masting")) +
  labs(y = "Oil content %", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(NA, 100))

pollData <- conifer[, c("mastEvent", "pollination")]
pollData$mastEvent <- as.factor(pollData$mastEvent)
pollData$pollination <- as.factor(pollData$pollination)
# Create ggplot object
poll <- ggplot(pollData, aes(x = pollination, fill = mastEvent)) +
  geom_bar(position = "dodge", colour = "black")  +
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#387E46", "1"="#C43142"), 
    labels = c("0"="No", "1"="Yes") 
  ) + scale_x_discrete(labels = c("animals" = "Animal", "wind" = "Wind", "wind and animals"="Both")) +
  labs(y = "Count", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(0, 35))

repData <- conifer[, c("mastEvent", "typeMonoOrDio")]
repData$mastEvent <- as.factor(repData$mastEvent)
repData$typeMonoOrDio <- as.factor(repData$typeMonoOrDio)
# Create ggplot object
rep <- ggplot(repData, aes(x = typeMonoOrDio, fill = mastEvent)) +
  geom_bar(position = "dodge", colour = "black") +
  annotate("text", x = 1.5, y = 55, 
           label = getAnnotation("Reproductive type", conifer_results),
           size = 2) +
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#387E46", "1"="#C43142"), 
    labels = c("0"="No", "1"="Yes") 
  ) +
  labs(y = "Count", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(0, 60))

droughtData <- conifer[, c("mastEvent", "droughtTolerance")]
droughtData$mastEvent <- as.factor(droughtData$mastEvent)
droughtData$droughtTolerance <- as.factor(droughtData$droughtTolerance)
# Create ggplot object
drought <- ggplot(droughtData, aes(x = droughtTolerance, fill = mastEvent)) +
  geom_bar(position = "dodge", colour = "black")  +
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#387E46", "1"="#C43142"), 
    labels = c("0"="No", "1"="Yes") 
  ) +
  labs(y = "Count", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(0, 25))

leafData <- conifer[, c("mastEvent", "leafLongevity")]
leafData$mastEvent <- as.factor(leafData$mastEvent)

leaf <- ggplot(leafData, aes(x = mastEvent, y = leafLongevity, fill = mastEvent)) +
  geom_violin(trim = FALSE, alpha = 0.5, colour = "black") +
  geom_point(position = position_jitter(width = 0.08), size = 1.5, alpha = 0.7)  +
  annotate("text", x = 1.5, y = 11, 
           label = getAnnotation("Leaf longevity (years)", conifer_results),
           size = 2) +
  scale_fill_manual(
    name = "Masting", 
    values = c("0"="#387E46", "1"="#C43142"), 
    labels = c("0"="No", "1"="Yes"))  + scale_x_discrete(labels = c("0" = "Non-masting", "1" = "Masting")) +
  labs(y = "Leaf longevity (year)", x = NULL) +
  theme_classic(base_size = 12) + scale_y_continuous(limits = c(NA, 12))

predation <- weight + fruit + oil + disp + dorm
pollination <- poll + rep + plot_layout(ncol = 3, widths = c(1,1,1))
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
