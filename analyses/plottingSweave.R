# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

setwd("C:/PhD/Project/PhD_thesis/mast_trait/")

source("analyses/dataCleaning.R")

library(ggplot2)
library(gridExtra)
library(grid)

custom_theme <- theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0.25),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.major = element_line(color = "gray85"),
    panel.grid.minor = element_blank()
  )

# remove rows with mastEvent data being NA
d <- d[!is.na(d$mastEvent), ]
# Make a subset of conifers only
conifer <- subset(d, familyName %in% c("Pinaceae","Taxodiaceae"))
# Make a subset of angiosperm only
angio <- subset(d, !(familyName %in% c("Pinaceae","Taxodiaceae")))

# Prepare a custom theme without legend
custom_theme_noleg <- custom_theme + theme(legend.position = "none")
# Extract legend from p_all
get_legend <- function(myplot) {
  g <- ggplotGrob(myplot)
  leg <- g$grobs[which(sapply(g$grobs, function(x) x$name) == "guide-box")]
  leg[[1]]
}


#seedDispersal
angio$seedDispersal <-  factor(angio$seedDispersal, levels = c("abiotic", "biotic","both"))
conifer$seedDispersal <-  factor(conifer$seedDispersal, levels = c("abiotic", "biotic","both"))
d$seedDispersal <-  factor(d$seedDispersal, levels = c("abiotic", "biotic","both"))

# Angiosperms plot
p_angio <- ggplot(angio, aes(x = mastEvent, fill = seedDispersal)) +
  geom_bar(position = "fill", color = "black", width = 0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("abiotic" = "#95B958", 
                               "biotic" = "#F4D166",
                               "both" = "#F59A3A")) +
  labs(title = "Angiosperms (N = 98)", x = "Masting", y = "Proportion") +
  custom_theme_noleg

# Conifers plot
p_conifer <- ggplot(conifer, aes(x = mastEvent, fill = seedDispersal)) +
  geom_bar(position = "fill", color = "black", width = 0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("abiotic" = "#95B958", 
                               "biotic" = "#F4D166",
                               "both" = "#F59A3A")) +
  labs(title = "Conifers (N = 61)", x = "Masting", y = "Proportion") +
  custom_theme_noleg

# All species plot (keep legend)
p_all <- ggplot(d, aes(x = mastEvent, fill = seedDispersal)) +
  geom_bar(position = "fill", color = "black", width = 0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("abiotic" = "#95B958", 
                               "biotic" = "#F4D166",
                               "both" = "#F59A3A")) +
  labs(title = "All Species (N = 159)", x = "Masting", y = "Proportion", fill = "Dispersal Mode") +
  custom_theme


shared_legend <- get_legend(p_all)

# Remove legend from main plots
p_all_noleg <- p_all + theme(legend.position = "none")

# Arrange plots in a row with legend on the right
combined <- arrangeGrob(
  p_angio, p_conifer, p_all_noleg,
  nrow = 1,
  widths = c(1,1,1)  # adjust widths if needed
)

#save to a pdf
pdf("output/figures/dispersalmode.pdf", width = 10, height = 4)
# Draw the combined plots and add the shared legend
grid.newpage()
grid.draw(arrangeGrob(combined, shared_legend, ncol=2, widths=c(3, 0.4)))
dev.off()

## seed weights
angioanimal <- angio[!is.na(angio$seedPredator) & angio$seedPredator != "no", ]
coniferanimal <- conifer[!is.na(conifer$seedPredator) & conifer$seedPredator != "no", ]
danimal <- d[!is.na(d$seedPredator) & d$seedPredator != "no", ]

angioanimal$logSeedWeights <- log10(angioanimal$seedWeights)

w_angio <- ggplot(angioanimal, aes(x = logSeedWeights, fill = mastEvent)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Y" = "#95B958", "N" = "#F4D166")) + labs(title = "Angiosperms (N = 98)", x = "Seed Weight (log10)", y = "Density") +
  custom_theme_noleg

coniferanimal$logSeedWeights <- log10(coniferanimal$seedWeights)

w_conifer <- ggplot(coniferanimal, aes(x = logSeedWeights, fill = mastEvent)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Y" = "#95B958", "N" = "#F4D166")) + labs(title = "Conifers (N = 61)", x = "Seed Weight (log10)", y = "Density") +
  custom_theme_noleg

danimal$logSeedWeights <- log10(danimal$seedWeights)

w_all <- ggplot(danimal, aes(x = logSeedWeights, fill = mastEvent)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Y" = "#95B958", "N" = "#F4D166")) + labs(title = "All Species (N = 159)", x = "Seed Weight (log10)", y = "Density") +
  custom_theme

shared_legend <- get_legend(w_all)

# Remove legend from main plots
w_all_noleg <- w_all + theme(legend.position = "none")

# Arrange plots in a row with legend on the right
combined <- arrangeGrob(
  w_angio, w_conifer, w_all_noleg,
  nrow = 1,
  widths = c(1,1,1)  # adjust widths if needed
)
#save to a pdf
pdf("output/figures/seedWeights.pdf", width = 10, height = 4)
# Draw the combined plots and add the shared legend
grid.newpage()
grid.draw(arrangeGrob(combined, shared_legend, ncol=2, widths=c(3, 0.4)))
dev.off()

# seed size

angioanimal$logseedSizeAve <- log10(angioanimal$seedSizeAve)

s_angio <- ggplot(angioanimal, aes(x = logseedSizeAve, fill = mastEvent)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Y" = "#95B958", "N" = "#F4D166")) + labs(title = "Angiosperms (N = 98)", x = "Seed Size (log10)", y = "Density") +
  custom_theme_noleg

coniferanimal$logseedSizeAve <- log10(coniferanimal$seedSizeAve)

s_conifer <- ggplot(coniferanimal, aes(x = logseedSizeAve, fill = mastEvent)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Y" = "#95B958", "N" = "#F4D166")) + labs(title = "Conifers (N = 61)", x = "Seed Size (log10)", y = "Density") +
  custom_theme_noleg

danimal$logseedSizeAve <- log10(danimal$seedSizeAve)

s_all <- ggplot(danimal, aes(x = logseedSizeAve, fill = mastEvent)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Y" = "#95B958", "N" = "#F4D166")) + labs(title = "All Species (N = 159)", x = "Seed Size (log10)", y = "Density") +
  custom_theme

shared_legend <- get_legend(s_all)

# Remove legend from main plots
s_all_noleg <- s_all + theme(legend.position = "none")

# Arrange plots in a row with legend on the right
combined <- arrangeGrob(
  s_angio, s_conifer, s_all_noleg,
  nrow = 1,
  widths = c(1,1,1)  # adjust widths if needed
)
#save to a pdf
pdf("output/figures/seedSize.pdf", width = 10, height = 4)
# Draw the combined plots and add the shared legend
grid.newpage()
grid.draw(arrangeGrob(combined, shared_legend, ncol=2, widths=c(3, 0.4)))
dev.off()

# Fruit size
angioanimal$logFruitSizeAve <- log10(angioanimal$fruitSizeAve)

f_angio <- ggplot(angioanimal, aes(x = logFruitSizeAve, fill = mastEvent)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Y" = "#95B958", "N" = "#F4D166")) + labs(title = "Angiosperms (N = 98)", x = "Fruit Size (log10)", y = "Density") +
  custom_theme_noleg

coniferanimal$logfruitSizeAve <- log10(coniferanimal$fruitSizeAve)

f_conifer <- ggplot(coniferanimal, aes(x = logfruitSizeAve, fill = mastEvent)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Y" = "#95B958", "N" = "#F4D166")) + labs(title = "Conifers (N = 61)", x = "Fruit Size (log10)", y = "Density") +
  custom_theme_noleg

danimal$logFruitSizeAve <- log10(danimal$fruitSizeAve)

f_all <- ggplot(danimal, aes(x = logFruitSizeAve, fill = mastEvent)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Y" = "#95B958", "N" = "#F4D166")) + labs(title = "All Species (N = 159)", x = "Fruit Size (log10)", y = "Density") +
  custom_theme

shared_legend <- get_legend(f_all)

# Remove legend from main plots
f_all_noleg <- f_all + theme(legend.position = "none")

# Arrange plots in a row with legend on the right
combined <- arrangeGrob(
  f_angio, f_conifer, f_all_noleg,
  nrow = 1,
  widths = c(1,1,1)  # adjust widths if needed
)
#save to a pdf
pdf("output/figures/fruitSize.pdf", width = 10, height = 4)
# Draw the combined plots and add the shared legend
grid.newpage()
grid.draw(arrangeGrob(combined, shared_legend, ncol=2, widths=c(3, 0.4)))
dev.off()

# Seed dormancy
angioanimal$seedDormancy <-  factor(angioanimal$seedDormancy, levels = c("Y", "N"))
coniferanimal$seedDormancy <-  factor(coniferanimal$seedDormancy, levels = c("Y", "N"))
danimal$seedDormancy <-  factor(danimal$seedDormancy, levels = c("Y", "N"))
# Angiosperms plot
d_angio <- ggplot(angioanimal, aes(x = mastEvent, fill = seedDormancy)) +
  geom_bar(position = "fill", color = "black", width = 0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("Y" = "#95B958", 
                               "N" = "#F4D166")) +
  labs(title = "Angiosperms (N = 98)", x = "Masting", y = "Proportion") +
  custom_theme_noleg

# Conifers plot
d_conifer <- ggplot(coniferanimal, aes(x = mastEvent, fill = seedDormancy)) +
  geom_bar(position = "fill", color = "black", width = 0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("Y" = "#95B958", 
                               "N" = "#F4D166")) +
  labs(title = "Conifers (N = 61)", x = "Masting", y = "Proportion") +
  custom_theme_noleg

# All species plot (keep legend)
d_all <- ggplot(danimal, aes(x = mastEvent, fill = seedDormancy)) +
  geom_bar(position = "fill", color = "black", width = 0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("Y" = "#95B958", 
                               "N" = "#F4D166")) +
  labs(title = "All Species (N = 159)", x = "Masting", y = "Proportion", fill = "seedDormancy") +
  custom_theme

shared_legend <- get_legend(d_all)

# Remove legend from main plots
d_all_noleg <- d_all + theme(legend.position = "none")

# Arrange plots in a row with legend on the right
combined <- arrangeGrob(
  d_angio, d_conifer, d_all_noleg,
  nrow = 1,
  widths = c(1,1,1)  # adjust widths if needed
)

#save to a pdf
pdf("output/figures/seedDormancy.pdf", width = 10, height = 4)
# Draw the combined plots and add the shared legend
grid.newpage()
grid.draw(arrangeGrob(combined, shared_legend, ncol=2, widths=c(3, 0.4)))
dev.off()


# pollination

# Reproductive mode
angio$typeMonoOrDio <- tolower(angio$typeMonoOrDio)
angio$typeMonoOrDio <- ifelse(angio$typeMonoOrDio == "monoecious", "Monoecious",
                              ifelse(angio$typeMonoOrDio == "dioecious", "Dioecious", 
                                     ifelse(angio$typeMonoOrDio == "polygamous", "Polygamous",NA)))
angio$typeMonoOrDio <- factor(angio$typeMonoOrDio, levels = c("Monoecious", "Dioecious","Polygamous"))
conifer$typeMonoOrDio <- tolower(conifer$typeMonoOrDio)
conifer$typeMonoOrDio <- ifelse(conifer$typeMonoOrDio == "monoecious", "Monoecious",
                              ifelse(conifer$typeMonoOrDio == "dioecious", "Dioecious", 
                                     ifelse(conifer$typeMonoOrDio == "polygamous", "Polygamous",NA)))
d$typeMonoOrDio <- factor(d$typeMonoOrDio, levels = c("Monoecious", "Dioecious","Polygamous"))
d$typeMonoOrDio <- tolower(d$typeMonoOrDio)
d$typeMonoOrDio <- ifelse(d$typeMonoOrDio == "monoecious", "Monoecious",
                                ifelse(d$typeMonoOrDio == "dioecious", "Dioecious", 
                                       ifelse(d$typeMonoOrDio == "polygamous", "Polygamous",NA)))
d$typeMonoOrDio <- factor(d$typeMonoOrDio, levels = c("Monoecious", "Dioecious","Polygamous"))
# Angiosperms plot
r_angio <- ggplot(angio, aes(x = mastEvent, fill = typeMonoOrDio)) +
  geom_bar(position = "fill", color = "black", width = 0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("Monoecious" = "#95B958", 
                               "Dioecious" = "#F4D166",
                               "Polygamous" = "#F59A3A")) +
  labs(title = "Angiosperms (N = 98)", x = "Masting", y = "Proportion") +
  custom_theme_noleg

# Conifers plot
r_conifer <- ggplot(conifer, aes(x = mastEvent, fill = typeMonoOrDio)) +
  geom_bar(position = "fill", color = "black", width = 0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("Monoecious" = "#95B958", 
                               "Dioecious" = "#F4D166",
                               "Polygamous" = "#F59A3A")) +
  labs(title = "Conifers (N = 61)", x = "Masting", y = "Proportion") +
  custom_theme_noleg

# All species plot (keep legend)
r_all <- ggplot(d, aes(x = mastEvent, fill = typeMonoOrDio)) +
  geom_bar(position = "fill", color = "black", width = 0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("Monoecious" = "#95B958", 
                               "Dioecious" = "#F4D166",
                               "Polygamous" = "#F59A3A")) +
  labs(title = "All Species (N = 159)", x = "Masting", y = "Proportion", fill = "typeMonoOrDio") +
  custom_theme

shared_legend <- get_legend(r_all)

# Remove legend from main plots
r_all_noleg <- r_all + theme(legend.position = "none")

# Arrange plots in a row with legend on the right
combined <- arrangeGrob(
  r_angio, r_conifer, r_all_noleg,
  nrow = 1,
  widths = c(1,1,1)  # adjust widths if needed
)
#save to a pdf
pdf("output/figures/reproductiveType.pdf", width = 10, height = 4)
# Draw the combined plots and add the shared legend
grid.newpage()
grid.draw(arrangeGrob(combined, shared_legend, ncol=2, widths=c(3, 0.4)))
dev.off()


# Drought tolerance
angio$droughtTolerance <- factor(angio$droughtTolerance, levels = c("High", "Moderate","Low"))
conifer$droughtTolerance <- factor(conifer$droughtTolerance, levels = c("High", "Moderate","Low"))
d$droughtTolerance <- factor(d$droughtTolerance, levels = c("High", "Moderate","Low"))
# Angiosperms plot
t_angio <- ggplot(angio, aes(x = mastEvent, fill = droughtTolerance)) +
  geom_bar(position = "fill", color = "black", width = 0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("High" = "#BA4C23", 
                               "Moderate" = "#F59A3A",
                               "Low" = "#FDC17B")) +
  labs(title = "Angiosperms (N = 98)", x = "Masting", y = "Proportion") +
  custom_theme_noleg

# Conifers plot
t_conifer <- ggplot(conifer, aes(x = mastEvent, fill = droughtTolerance)) +
  geom_bar(position = "fill", color = "black", width = 0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("High" = "#BA4C23", 
                               "Moderate" = "#F59A3A",
                               "Low" = "#FDC17B")) +
  labs(title = "Conifers (N = 61)", x = "Masting", y = "Proportion") +
  custom_theme_noleg

# All species plot (keep legend)
t_all <- ggplot(d, aes(x = mastEvent, fill = droughtTolerance)) +
  geom_bar(position = "fill", color = "black", width = 0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("High" = "#BA4C23", 
                               "Moderate" = "#F59A3A",
                               "Low" = "#FDC17B")) +
  labs(title = "All Species (N = 159)", x = "Masting", y = "Proportion", fill = "droughtTolerance") +
  custom_theme

shared_legend <- get_legend(t_all)

# Remove legend from main plots
t_all_noleg <- t_all + theme(legend.position = "none")

# Arrange plots in a row with legend on the right
combined <- arrangeGrob(
  t_angio, t_conifer, t_all_noleg,
  nrow = 1,
  widths = c(1,1,1)  # adjust widths if needed
)

#save to a pdf
pdf("output/figures/droughtTolerance.pdf", width = 10, height = 4)
# Draw the combined plots and add the shared legend
grid.newpage()
grid.draw(arrangeGrob(combined, shared_legend, ncol=2, widths=c(3, 0.4)))
dev.off()

# flowering duration
angio$floweringDuration <- factor(angio$floweringDuration, levels = c("<10", "10-30",">30"))
conifer$floweringDuration <- factor(conifer$floweringDuration, levels = c("<10", "10-30",">30"))
d$floweringDuration <- factor(d$floweringDuration, levels = c("<10", "10-30",">30"))
# Angiosperms plot
f_angio <- ggplot(angio, aes(x = mastEvent, fill = floweringDuration)) +
  geom_bar(position = "fill", color = "black", width = 0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c(">30" = "#BA4C23", 
                               "10-30" = "#F59A3A",
                               "<10" = "#FDC17B")) +
  labs(title = "Angiosperms (N = 98)", x = "Masting", y = "Proportion") +
  custom_theme_noleg

# Conifers plot
f_conifer <- ggplot(conifer, aes(x = mastEvent, fill = floweringDuration)) +
  geom_bar(position = "fill", color = "black", width = 0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c(">30" = "#BA4C23", 
                               "10-30" = "#F59A3A",
                               "<10" = "#FDC17B")) +
  labs(title = "Conifers (N = 61)", x = "Masting", y = "Proportion") +
  custom_theme_noleg

# All species plot (keep legend)
f_all <- ggplot(d, aes(x = mastEvent, fill = floweringDuration)) +
  geom_bar(position = "fill", color = "black", width = 0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c(">30" = "#BA4C23", 
                               "10-30" = "#F59A3A",
                               "<10" = "#FDC17B")) +
  labs(title = "All Species (N = 159)", x = "Masting", y = "Proportion", fill = "floweringDuration") +
  custom_theme

shared_legend <- get_legend(f_all)

# Remove legend from main plots
f_all_noleg <- f_all + theme(legend.position = "none")

# Arrange plots in a row with legend on the right
combined <- arrangeGrob(
  f_angio, f_conifer, f_all_noleg,
  nrow = 1,
  widths = c(1,1,1)  # adjust widths if needed
)

#save to a pdf
pdf("output/figures/flowerDuration.pdf", width = 10, height = 4)
# Draw the combined plots and add the shared legend
grid.newpage()
grid.draw(arrangeGrob(combined, shared_legend, ncol=2, widths=c(3, 0.4)))
dev.off()

# Oil content
o_angio <- ggplot(angioanimal, aes(x = oilContent, fill = mastEvent)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Y" = "#95B958", "N" = "#F4D166")) + labs(title = "Angiosperms (N = 98)", x = "Oil Content", y = "Density") +
  custom_theme_noleg

o_conifer <- ggplot(coniferanimal, aes(x = oilContent, fill = mastEvent)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Y" = "#95B958", "N" = "#F4D166")) + labs(title = "Conifers (N = 61)", x = "Oil Content", y = "Density") +
  custom_theme_noleg

o_all <- ggplot(danimal, aes(x = oilContent, fill = mastEvent)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Y" = "#95B958", "N" = "#F4D166")) + labs(title = "All Species (N = 159)", x = "Oil Content", y = "Density") +
  custom_theme

shared_legend <- get_legend(o_all)

# Remove legend from main plots
o_all_noleg <- o_all + theme(legend.position = "none")

# Arrange plots in a row with legend on the right
combined <- arrangeGrob(
  o_angio, o_conifer, o_all_noleg,
  nrow = 1,
  widths = c(1,1,1)  # adjust widths if needed
)

#save to a pdf
pdf("output/figures/oilContent.pdf", width = 10, height = 4)
# Draw the combined plots and add the shared legend
grid.newpage()
grid.draw(arrangeGrob(combined, shared_legend, ncol=2, widths=c(3, 0.4)))
dev.off()

# Leaf longevity

l_angio <- ggplot(angio, aes(x = leafLongevity, fill = mastEvent)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Y" = "#95B958", "N" = "#F4D166")) + labs(title = "Angiosperms (N = 98)", x = "Leaf Longevity", y = "Density") +
  custom_theme_noleg

l_conifer <- ggplot(conifer, aes(x = leafLongevity, fill = mastEvent)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Y" = "#95B958", "N" = "#F4D166")) + labs(title = "Conifers (N = 61)", x = "Leaf Longevity", y = "Density")  +
  custom_theme

shared_legend <- get_legend(l_conifer)

# Remove legend from main plots
l_conifer_noleg <- l_conifer + theme(legend.position = "none")

# Arrange plots in a row with legend on the right
combined <- arrangeGrob(
  l_angio, l_conifer_noleg,
  nrow = 1,
  widths = c(1,1)  # adjust widths if needed
)

#save to a pdf
pdf("output/figures/leafLongevity.pdf", width = 10, height = 4)
# Draw the combined plots and add the shared legend
grid.newpage()
grid.draw(arrangeGrob(combined, shared_legend, ncol=2, widths=c(3, 0.4)))
dev.off()

# pollination
angio$pollination <- factor(angio$pollination, levels = c("wind", "animals", "wind and animals"))
p_angio <- ggplot(angio, aes(x = mastEvent, fill = pollination)) +
  geom_bar(position = "fill", color = "black", width=0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("wind" = "#95B958", "animals" = "#F4D166", "wind and animals" = "#6194BF")) +
  labs(
    title = "Angiosperms (N = 98)",
    x = "Masting",
    y = "Proportion"  ) +
  custom_theme_noleg
conifer$pollination <- factor(conifer$pollination, levels = c("wind", "animals", "wind and animals"))
p_conifer <- ggplot(conifer, aes(x = mastEvent, fill = pollination)) +
  geom_bar(position = "fill", color = "black", width=0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("wind" = "#95B958", "animals" = "#F4D166", "wind and animals" = "#6194BF")) +
  labs(
    title = "Conifers (N = 61)",
    x = "Masting",
    y = "Proportion"  ) +
  custom_theme_noleg
d$pollination <- factor(d$pollination, levels = c("wind", "animals", "wind and animals"))
p_all <- ggplot(d, aes(x = mastEvent, fill = pollination)) +
  geom_bar(position = "fill", color = "black", width=0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("wind" = "#95B958", "animals" = "#F4D166", "wind and animals" = "#6194BF")) +
  labs(
    title = "All Species (N = 159)",
    x = "Masting",
    y = "Proportion", fill ="pollination"  ) +
  custom_theme

shared_legend <- get_legend(p_all)

# Remove legend from main plots
p_all_noleg <- p_all + theme(legend.position = "none")

# Arrange plots in a row with legend on the right
combined <- arrangeGrob(
  p_angio, p_conifer, p_all_noleg,
  nrow = 1,
  widths = c(1,1,1)  # adjust widths if needed
)

#save to a pdf
pdf("output/figures/pollination.pdf", width = 10, height = 4)
# Draw the combined plots and add the shared legend
grid.newpage()
grid.draw(arrangeGrob(combined, shared_legend, ncol=2, widths=c(3, 0.4)))
dev.off()