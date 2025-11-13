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

# Prepare a custom theme without legend
custom_theme_noleg <- custom_theme + theme(legend.position = "none")

# Angiosperms plot
p_angio <- ggplot(angio, aes(x = mastEvent, fill = seedDispersal)) +
  geom_bar(position = "fill", color = "black", width = 0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("abiotic" = "#95B958", 
                               "biotic" = "#F4D166",
                               "both" = "lightblue")) +
  labs(title = "Angiosperms", x = "Masting", y = "Proportion") +
  custom_theme_noleg

# Conifers plot
p_conifer <- ggplot(conifer, aes(x = mastEvent, fill = seedDispersal)) +
  geom_bar(position = "fill", color = "black", width = 0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("abiotic" = "#95B958", 
                               "biotic" = "#F4D166",
                               "both" = "lightblue")) +
  labs(title = "Conifers", x = "Masting", y = "Proportion") +
  custom_theme_noleg

# All species plot (keep legend)
p_all <- ggplot(d, aes(x = mastEvent, fill = seedDispersal)) +
  geom_bar(position = "fill", color = "black", width = 0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("abiotic" = "#95B958", 
                               "biotic" = "#F4D166",
                               "both" = "lightblue")) +
  labs(title = "All Species", x = "Masting", y = "Proportion", fill = "Dispersal Mode") +
  custom_theme

# Extract legend from p_all
get_legend <- function(myplot) {
  g <- ggplotGrob(myplot)
  leg <- g$grobs[which(sapply(g$grobs, function(x) x$name) == "guide-box")]
  leg[[1]]
}
shared_legend <- get_legend(p_all)

# Remove legend from main plots
p_all_noleg <- p_all + theme(legend.position = "none")

# Arrange plots in a row with legend on the right
combined <- arrangeGrob(
  p_angio, p_conifer, p_all_noleg,
  nrow = 1,
  widths = c(1,1,1)  # adjust widths if needed
)

# Draw the combined plots and add the shared legend
grid.newpage()
grid.draw(arrangeGrob(combined, shared_legend, ncol=2, widths=c(3, 0.4)))

## seed weights
angiobiotic <- angio[angio$seedDispersal == c("biotic","both"), ]
coniferbiotic <- conifer[conifer$seedDispersal == c("biotic","both"), ]
dbiotic <- d[d$seedDispersal == c("biotic","both"), ]
angio$logSeedWeights <- log10(angio$seedWeights)

w_angio <- ggplot(angio, aes(x = logSeedWeights, fill = mastEvent)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Y" = "#95B958", "N" = "#F4D166")) + labs(title = "Angiosperms", x = "Seed Weight (log10)", y = "Density") +
  custom_theme_noleg

conifer$logSeedWeights <- log10(conifer$seedWeights)

w_conifer <- ggplot(conifer, aes(x = logSeedWeights, fill = mastEvent)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Y" = "#95B958", "N" = "#F4D166")) + labs(title = "Conifers", x = "Seed Weight (log10)", y = "Density") +
  custom_theme_noleg

d$logSeedWeights <- log10(d$seedWeights)

w_all <- ggplot(d, aes(x = logSeedWeights, fill = mastEvent)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Y" = "#95B958", "N" = "#F4D166")) + labs(title = "AllSpecies", x = "Seed Weight (log10)", y = "Density") +
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

# Draw the combined plots and add the shared legend
grid.newpage()
grid.draw(arrangeGrob(combined, shared_legend, ncol=2, widths=c(3, 0.4)))

# seed size
