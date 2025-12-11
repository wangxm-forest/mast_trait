#### Phylogeny logistic regression ####
## Started by Mao ##
## Nov-28-2025 ##

library(caper)
library(geiger)
library(phylolm)
library(brglm2)
library(detectseparation)

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
library(brglm2)

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
