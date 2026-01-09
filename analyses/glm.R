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

rm(list = ls())
options(stringsAsFactors = FALSE)

setwd("C:/PhD/Project/PhD_thesis/mast_trait")


# Load & prepare Data
d <- read.csv("data/cleanSilvics.csv")
phytree <- read.tree("output/silvicsPhylogenyFull.tre")

d$latbi <- gsub(" ", "_", d$latbi)
d <- d[!is.na(d$mastEvent), ]
d$mastEvent <- ifelse(d$mastEvent == "Y", 1, 0)

conifer <- d[d$familyName %in% c("Pinaceae", "Taxodiaceae"), ]
angio   <- d[!(d$familyName %in% c("Pinaceae", "Taxodiaceae")), ]

# log10-transform + scale
log_scale <- function(x) scale(log10(x))

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

# Continuous traits
conifer$logSeedWeightStd <- log_scale(conifer$seedWeights)
angio$logSeedWeightStd   <- log_scale(angio$seedWeights)

conifer$logFruitStd <- log_scale(conifer$fruitSizeAve)
angio$logFruitStd   <- log_scale(angio$fruitSizeAve)

conifer$logSeedSizeStd <- log_scale(conifer$seedSizeAve)
angio$logSeedSizeStd   <- log_scale(angio$seedSizeAve)

phytree$node.label <- NULL

phyconifer <- drop.tip(phytree, setdiff(phytree$tip.label, conifer$latbi))
phyangio   <- drop.tip(phytree, setdiff(phytree$tip.label, angio$latbi))

rownames(conifer) <- conifer$latbi
rownames(angio)   <- angio$latbi



# Make some functions

# Extract key results from a phyloglm object
tidy_phyloglm <- function(model) {
  s <- summary(model)
  coefs <- s$coefficients
  out <- data.frame(
    term = rownames(coefs),
    estimate = coefs[, 1],
    std_error = coefs[, 2],
    z_value = coefs[, 3],
    p_value = coefs[, 4],
    alpha = model$alpha,
    row.names = NULL
  )
  return(out)
}


# Run model and attach metadata
run_model <- function(formula, data, phy, method, trait_name) {
  model <- phyloglm(formula, phy = phy, data = data, method = method)
  tbl <- tidy_phyloglm(model)
  tbl$trait <- trait_name
  tbl$model_formula <- deparse(formula)
  return(tbl)
}


# Prepare Trait Columns
categorical_traits <- c("seedDispersal", "seedDormancy",
                        "typeMonoOrDio", "droughtTolerance", "pollination")

for (t in categorical_traits) {
  conifer[[t]] <- factor(conifer[[t]])
  angio[[t]]   <- factor(angio[[t]])
}




# Model definitions
model_list <- list(
  list(name="Seed dispersal (conifer)", formula=mastEvent ~ seedDispersal, data=conifer, phy=phyconifer, method="logistic_MPLE"),
  list(name="Seed dispersal (angio)",   formula=mastEvent ~ seedDispersal, data=angio,   phy=phyangio,   method="logistic_MPLE"),
  
  list(name="Pollination (angio)",   formula=mastEvent ~ pollination, data=angio,   phy=phyangio,   method="logistic_MPLE"),
  
  list(name="Seed dormancy (conifer)",  formula=mastEvent ~ seedDormancy, data=conifer, phy=phyconifer, method="logistic_MPLE"),
  list(name="Seed dormancy (angio)",    formula=mastEvent ~ seedDormancy, data=angio,   phy=phyangio,   method="logistic_MPLE"),
  
  list(name="Mono/Dio (conifer)",       formula=mastEvent ~ typeMonoOrDio, data=conifer, phy=phyconifer, method="logistic_MPLE"),
  list(name="Mono/Dio (angio)",         formula=mastEvent ~ typeMonoOrDio, data=angio,   phy=phyangio,   method="logistic_MPLE"),
  
  list(name="Seed weight (conifer)",    formula=mastEvent ~ logSeedWeightStd, data=conifer, phy=phyconifer, method="logistic_IG10"),
  list(name="Seed weight (angio)",      formula=mastEvent ~ logSeedWeightStd, data=angio,   phy=phyangio,   method="logistic_IG10"),
  
  list(name="Fruit size (conifer)",     formula=mastEvent ~ logFruitStd, data=conifer, phy=phyconifer, method="logistic_IG10"),
  list(name="Fruit size (angio)",       formula=mastEvent ~ logFruitStd, data=angio,   phy=phyangio,   method="logistic_IG10"),
  
  list(name="Seed size (conifer)",      formula=mastEvent ~ logSeedSizeStd, data=conifer, phy=phyconifer, method="logistic_IG10"),
  list(name="Seed size (angio)",        formula=mastEvent ~ logSeedSizeStd, data=angio,   phy=phyangio,   method="logistic_IG10"),
  
  list(name="Oil content (angio)",      formula=mastEvent ~ oilContent, data=angio, phy=phyangio, method="logistic_IG10"),
  
  list(name="Leaf longevity (conifer)", formula=mastEvent ~ leafLongevity, data=conifer, phy=phyconifer, method="logistic_IG10"),
  list(name="Leaf longevity (angio)",   formula=mastEvent ~ leafLongevity, data=angio,   phy=phyangio,   method="logistic_IG10"),
  
  list(name="Drought tol (conifer)",    formula=mastEvent ~ droughtTolerance, data=conifer, phy=phyconifer, method="logistic_IG10"),
  list(name="Drought tol (angio)",      formula=mastEvent ~ droughtTolerance, data=angio,   phy=phyangio,   method="logistic_IG10")
)


# Run models

results <- NULL

for (m in model_list) {
  cat("Running model:", m$name, "\n")
  
  tbl <- run_model(
    formula = m$formula,
    data    = m$data,
    phy     = m$phy,
    method  = m$method,
    trait_name = m$name
  )
  
  results <- rbind(results, tbl)
}


#write.csv(results, "output/pglsResults.csv", row.names = FALSE)

#Make a table to present the results:
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
  "z_value",
  "p_value",
  "signif",
  "alpha"
)]

colnames(final_table) <- c(
  "Trait (group)",
  "Predictor",
  "Estimate",
  "SE",
  "Z",
  "P",
  "Sig.",
  "Phylo α"
)


final_table$Predictor <- gsub("logFruitStd",     "Fruit size (log, std)", final_table$Predictor)
final_table$Predictor <- gsub("logSeedWeightStd","Seed weight (log, std)", final_table$Predictor)
final_table$Predictor <- gsub("logSeedSizeStd",  "Seed size (log, std)", final_table$Predictor)
final_table$Predictor <- gsub(,  "Biotic dispersed (compared to Abiotic)", final_table$Predictor)
final_table$Predictor <- gsub("seedDispersalboth",  "Abiotic and Biotic dispersed (compared to Abiotic)", final_table$Predictor)
final_table$Predictor <- gsub("pollinationwind",  "Wind pollinated (compared to Animal pollinated)", final_table$Predictor)
final_table$Predictor <- gsub("pollinationwind and animals (compared to Animal pollinated)",  "Animal and wind pollinated", final_table$Predictor)
final_table$Predictor <- gsub("seedDormancyY",  "Dormant", final_table$Predictor)
final_table$Predictor <- gsub("typeMonoOrDioMonoecious",  "Monoecious (compared to Dioecious)", final_table$Predictor)
final_table$Predictor <- gsub("typeMonoOrDioPolygamous",  "Polygamous (compared to Dioecious)", final_table$Predictor)
final_table$Predictor <- gsub("oilContent",  "Oil content", final_table$Predictor)
final_table$Predictor <- gsub("leafLongevity",  "Leaf longevity", final_table$Predictor)
final_table$Predictor <- gsub("droughtToleranceLow",  "Low drought tolerated (compared to High drought tolerated)", final_table$Predictor)
final_table$Predictor <- gsub("droughtToleranceModerate",  "Moderate drought tolerated (compared to High drought tolerated)", final_table$Predictor)


rownames(final_table) <- NULL
table_grob <- tableGrob(
  final_table,
  rows = NULL,
  theme = ttheme_minimal(
    core = list(fg_params = list(fontsize = 9)),
    colhead = list(fg_params = list(fontsize = 10, fontface = "bold"))
  )
)
ggsave(
  filename = "output/phyloglmResultsTable.pdf",
  plot = table_grob,
  width = 10,
  height = 8
)

####Use conifer/angio as a fixed effect in the model####

# Create group variable
d$group <- ifelse(d$familyName %in% c("Pinaceae", "Taxodiaceae"),
                  "conifer", "angiosperm")
d$group <- factor(d$group)

# Prepare Trait Columns
categorical_traits <- c("seedDispersal", "seedDormancy",
                        "typeMonoOrDio", "droughtTolerance", "pollination")

for (t in categorical_traits) {
  d[[t]] <- factor(d[[t]])
}

# Continuous traits
d$logSeedWeightStd <- log_scale(d$seedWeights)
d$logFruitStd <- log_scale(d$fruitSizeAve)
d$logSeedSizeStd   <- log_scale(d$seedSizeAve)

# Extract key results from a phyloglm object
tidy_glm <- function(model) {
  s <- summary(model)$coefficients
  out <- data.frame(
    term = rownames(s),
    estimate = s[,1],
    std_error = s[,2],
    z_value = s[,3],
    p_value = s[,4],
    row.names = NULL
  )
  return(out)
}

# Run model and attach metadata
run_model <- function(formula, data, method, trait_name) {
  model <- glm(formula, data = data, method = method)
  tbl <- tidy_glm(model)
  tbl$trait <- trait_name
  tbl$model_formula <- deparse(formula)
  return(tbl)
}

# Model definitions
model_list <- list(
  list(name="Seed dispersal", formula=mastEvent ~ seedDispersal + group, data=d, method="brglmFit"),
 
  list(name="Pollination", formula=mastEvent ~ pollination + group, data=d, method="brglmFit"),
  
  list(name="Seed dormancy", formula=mastEvent ~ seedDormancy + group, data=d, method="brglmFit"),
 
  list(name="Mono/Dio", formula=mastEvent ~ typeMonoOrDio + group, data=d, method="brglmFit"),
 
  list(name="Seed weight (log)",    formula=mastEvent ~ logSeedWeightStd + group, data=d, method="brglmFit"),
  
  list(name="Fruit size (log)",     formula=mastEvent ~ logFruitStd + group, data=d, method="brglmFit"),
  
  list(name="Seed size (log)",      formula=mastEvent ~ logSeedSizeStd + group, data=d, method="brglmFit"),
  
  list(name="Oil content %",      formula=mastEvent ~ oilContent + group, data=d, method="brglmFit"),
  
  list(name="Leaf longevity (years)", formula=mastEvent ~ leafLongevity + group, data=d, method="brglmFit"),
  
  list(name="Drought tolerance",    formula=mastEvent ~ droughtTolerance + group, data=d, method="brglmFit")
)

# Run models

results_glm <- NULL

for (m in model_list) {
  cat("Running model:", m$name, "\n")
  
  tbl <- run_model(
    formula = m$formula,
    data    = m$data,
    method  = m$method,
    trait_name = m$name
  )
  
  results_glm <- rbind(results_glm, tbl)
}

#write.csv(results_glm, "output/glmResults.csv", row.names = FALSE)

####Plot the probability for all traits####
invlogit <- function(x) 1 / (1 + exp(-x))

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

plot_prob_with_signif <- function(pred_df, trait_name, results_df, trait_var) {
  
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
  ggplot(pred_df, aes(x = trait_level, y = prob, color = group, group = group)) +
    geom_point(size = 3) +
    geom_text(aes(label = signif), vjust = -0.5, size = 5, show.legend = FALSE,
              position = position_dodge(width = 1)) +
    labs(
      x = trait_name,
      y = "Predicted probability of masting"
    ) + scale_color_manual(values = my_colors) +
    custom_theme + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

# Pollination
poll <- subset(results_glm, trait == "Pollination")
poll

coefs <- setNames(poll$estimate, poll$term)
coefs

predPoll <- expand.grid(
  trait_level = c("animal", "wind", "wind and animals"),
  group = c("angiosperm", "conifer")
)
predPoll

predPoll$logit <- coefs["(Intercept)"]

predPoll$logit <- predPoll$logit +
  ifelse(predPoll$group == "conifer", coefs["groupconifer"], 0)

predPoll$logit <- predPoll$logit +
  ifelse(predPoll$trait_level == "wind", coefs["pollinationwind"], 0) +
  ifelse(predPoll$trait_level == "wind and animals", coefs["pollinationwind and animals"], 0)

predPoll$prob <- invlogit(predPoll$logit)
predPoll

# Seed dispersal
disp <- subset(results_glm, trait == "Seed dispersal")
disp

coefs <- setNames(disp$estimate, disp$term)
coefs

predDisp <- expand.grid(
  trait_level = c("abiotic", "biotic", "both"),
  group = c("angiosperm", "conifer")
)
predDisp

predDisp$logit <- coefs["(Intercept)"]

predDisp$logit <- predDisp$logit +
  ifelse(predDisp$group == "conifer", coefs["groupconifer"], 0)

predDisp$logit <- predDisp$logit +
  ifelse(predDisp$trait_level == "biotic", coefs["seedDispersalbiotic"], 0) +
  ifelse(predDisp$trait_level == "both", coefs["seedDispersalboth"], 0)

predDisp$prob <- invlogit(predDisp$logit)
predDisp

# Seed dormancy
dorm <- subset(results_glm, trait == "Seed dormancy")
dorm

coefs <- setNames(dorm$estimate, dorm$term)
coefs

predDorm <- expand.grid(
  trait_level = c("N","Y"),
  group = c("angiosperm", "conifer")
)
predDorm

predDorm$logit <- coefs["(Intercept)"]

predDorm$logit <- predDorm$logit +
  ifelse(predDorm$group == "conifer", coefs["groupconifer"], 0)

predDorm$logit <- predDorm$logit +
  ifelse(predDorm$trait_level == "Y", coefs["seedDormancyY"], 0)

predDorm$prob <- invlogit(predDorm$logit)
predDorm

# Mono/Dio
mono <- subset(results_glm, trait == "Mono/Dio")
mono

coefs <- setNames(mono$estimate, mono$term)
coefs

predMono <- expand.grid(
  trait_level = c("Dioecious", "Monoecious", "Polygamous"),
  group = c("angiosperm", "conifer")
)
predMono

predMono$logit <- coefs["(Intercept)"]

predMono$logit <- predMono$logit +
  ifelse(predMono$group == "conifer", coefs["groupconifer"], 0)

predMono$logit <- predMono$logit +
  ifelse(predMono$trait_level == "Monoecious", coefs["typeMonoOrDioMonoecious"], 0)+
  ifelse(predMono$trait_level == "Polygamous", coefs["typeMonoOrDioPolygamous"], 0)

predMono$prob <- invlogit(predMono$logit)
predMono

# Drought tolerance
drought <- subset(results_glm, trait == "Drought tolerance")
drought

coefs <- setNames(drought$estimate, drought$term)
coefs

preddrought <- expand.grid(
  trait_level = c("High", "Low", "Moderate"),
  group = c("angiosperm", "conifer")
)
preddrought

preddrought$logit <- coefs["(Intercept)"]

preddrought$logit <- preddrought$logit +
  ifelse(preddrought$group == "conifer", coefs["groupconifer"], 0)

preddrought$logit <- preddrought$logit +
  ifelse(preddrought$trait_level == "Low", coefs["droughtToleranceLow"], 0)+
  ifelse(preddrought$trait_level == "Moderate", coefs["droughtToleranceModerate"], 0)

preddrought$prob <- invlogit(preddrought$logit)
preddrought

#Plotting
plotPoll <- plot_prob_with_signif(predPoll, 
                      trait_name = "Pollination", 
                      results_df = results_glm, 
                      trait_var = "pollination") + theme(legend.position = "none")

plotDisp <- plot_prob_with_signif(predDisp, 
                      trait_name = "Seed dispersal", 
                      results_df = results_glm, 
                      trait_var = "seedDispersal") + theme(legend.position = "none")

plotDorm <- plot_prob_with_signif(predDorm, 
                      trait_name = "Seed dormancy", 
                      results_df = results_glm, 
                      trait_var = "seedDormancy") + theme(legend.position = "none")

plotMono <- plot_prob_with_signif(predMono, 
                      trait_name = "Mono/Dio", 
                      results_df = results_glm, 
                      trait_var = "typeMonoOrDio") + theme(legend.position = "none")

plotDrought <- plot_prob_with_signif(predMono, 
                      trait_name = "Drought tolerance", 
                      results_df = results_glm, 
                      trait_var = "droughtTolerance")
shared_legend <- get_legend(plotDrought)
plotDrought <- plot_prob_with_signif(predMono, 
                                     trait_name = "Drought tolerance", 
                                     results_df = results_glm, 
                                     trait_var = "droughtTolerance") + theme(legend.position = "none")

pdf("output/figures/glmCat.pdf", width = 10, height = 10)
grid.arrange(plotPoll, plotDisp, plotDorm, plotMono,plotDrought, shared_legend, nrow = 2, ncol = 3)
dev.off()
## Continuous traits

cont_effects <- results_glm %>%
  filter(trait %in% c("Seed weight (log)", "Fruit size (log)", "Seed size (log)","Leaf longevity (years)", "Oil content %"))%>%
  filter(!term %in% "(Intercept)")

cont_effects <- cont_effects %>%
  mutate(
    lower = estimate - 1.96 * std_error,
    upper = estimate + 1.96 * std_error,
    signif = ifelse(p_value < 0.05, "*", "ns"),
    type = ifelse(grepl("group", term), "Group effect", "Trait effect"),
    group = ifelse(grepl("groupconifer", term), "conifer", "angiosperm")
  )
pdf("output/figures/glmCon.pdf", width = 10, height = 10)
ggplot(cont_effects, aes(x = trait, y = estimate, color = group)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2,
                position = position_dodge(width = 0.5)) +
  geom_text(aes(label = signif), vjust = -1, size = 5, show.legend = FALSE,
            position = position_dodge(width = 1)) +
  labs(
    x = "Trait / Group",
    y = "Effect size (log-odds)"
  ) +
  custom_theme + scale_color_manual(values = my_colors) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

# Plot the mean and SE
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
  "Trait", "Term", "Estimate", "Std_Error", "P_value", "Lambda", "N"
)
colnames(resultsAngio) <- c(
  "Trait", "Term", "Estimate", "Std_Error", "P_value", "Lambda", "N"
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
dev.off()

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
dev.off()
