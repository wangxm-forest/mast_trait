#### Phylogeny logistic regression ####
## Started by Mao ##
## Nov-28-2025 ##

library(caper)
library(geiger)
library(phylolm)
library(brglm2)
library(detectseparation)
library(gridExtra)
library(grid)

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

# log10-transform + scale
log_scale <- function(x) scale(log10(x))


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

# Continuous traits
conifer$logSeedWeightStd <- log_scale(conifer$seedWeights)
angio$logSeedWeightStd   <- log_scale(angio$seedWeights)

conifer$logFruitStd <- log_scale(conifer$fruitSizeAve)
angio$logFruitStd   <- log_scale(angio$fruitSizeAve)

conifer$logSeedSizeStd <- log_scale(conifer$seedSizeAve)
angio$logSeedSizeStd   <- log_scale(angio$seedSizeAve)


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


write.csv(results, "output/pglsResults.csv", row.names = FALSE)


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
 
  list(name="Seed weight",    formula=mastEvent ~ logSeedWeightStd + group, data=d, method="brglmFit"),
  
  list(name="Fruit size",     formula=mastEvent ~ logFruitStd + group, data=d, method="brglmFit"),
  
  list(name="Seed size",      formula=mastEvent ~ logSeedSizeStd + group, data=d, method="brglmFit"),
  
  list(name="Oil content",      formula=mastEvent ~ oilContent + group, data=d, method="brglmFit"),
  
  list(name="Leaf longevity", formula=mastEvent ~ leafLongevity + group, data=d, method="brglmFit"),
  
  list(name="Drought tolerance",    formula=mastEvent ~ droughtTolerance + group, data=d, method="brglmFit")
)

# Run models

results <- NULL

for (m in model_list) {
  cat("Running model:", m$name, "\n")
  
  tbl <- run_model(
    formula = m$formula,
    data    = m$data,
    method  = m$method,
    trait_name = m$name
  )
  
  results <- rbind(results, tbl)
}

write.csv(results, "output/glmResults.csv", row.names = FALSE)

####Plot the probability for all traits####
invlogit <- function(x) 1 / (1 + exp(-x))

my_colors <- c("angiosperm" = "#95B958",
               "conifer"    = "#6194BF")

custom_theme <- theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0.25),
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
poll <- subset(results, trait == "Pollination")
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
disp <- subset(results, trait == "Seed dispersal")
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
dorm <- subset(results, trait == "Seed dormancy")
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
mono <- subset(results, trait == "Mono/Dio")
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
drought <- subset(results, trait == "Drought tolerance")
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
                      results_df = results, 
                      trait_var = "pollination") + theme(legend.position = "none")

plotDisp <- plot_prob_with_signif(predDisp, 
                      trait_name = "Seed dispersal", 
                      results_df = results, 
                      trait_var = "seedDispersal") + theme(legend.position = "none")

plotDorm <- plot_prob_with_signif(predDorm, 
                      trait_name = "Seed dormancy", 
                      results_df = results, 
                      trait_var = "seedDormancy") + theme(legend.position = "none")

plotMono <- plot_prob_with_signif(predMono, 
                      trait_name = "Mono/Dio", 
                      results_df = results, 
                      trait_var = "typeMonoOrDio") + theme(legend.position = "none")

plotDrought <- plot_prob_with_signif(predMono, 
                      trait_name = "Drought tolerance", 
                      results_df = results, 
                      trait_var = "droughtTolerance")
shared_legend <- get_legend(plotDrought)
plotDrought <- plot_prob_with_signif(predMono, 
                                     trait_name = "Drought tolerance", 
                                     results_df = results, 
                                     trait_var = "droughtTolerance") + theme(legend.position = "none")

pdf("output/figures/glmCat.pdf", width = 10, height = 10)
grid.arrange(plotPoll, plotDisp, plotDorm, plotMono,plotDrought, shared_legend, nrow = 2, ncol = 3)
dev.off()
## Continuous traits

cont_effects <- results %>%
  filter(trait %in% c("Seed weight", "Fruit size", "Seed size","Leaf longevity"))%>%
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
