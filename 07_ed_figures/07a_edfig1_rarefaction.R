# ED Fig 1: Rarefaction analysis (observed vs rarefied species richness)
# 数据: SI Table 3 (每样本 Bacteria/Archaea/Eukaryotes 的测序深度、观测丰富度与稀释丰富度)
# 输出: outputs/07a_edfig1_observed_vs_rarefied.pdf/.png
# 运行: Rscript 07_ed_figures/07a_edfig1_rarefaction.R
library(openxlsx)
library(ggplot2)

t3 <- read.xlsx("data/Supplementary_Tables.xlsx", sheet = "Table 3",
                startRow = 3, check.names = FALSE)

plot_data <- data.frame(
  S = t3$Observed_richness,
  Srare = t3$Rarefied_richness
)

model <- lm(Srare ~ S, data = plot_data)
rsq <- summary(model)$r.squared

p <- ggplot(plot_data, aes(x = S, y = Srare)) +
  geom_point(color = "#3C6E9E", size = 0.5) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "black",
    linewidth = 0.5
  ) +
  labs(
    x = "Observed No. of Species",
    y = "Rarefied No. of Species"
  ) +
  annotate(
    "text",
    x = max(plot_data$S),
    y = max(plot_data$Srare),
    label = paste0("R² = ", round(rsq, 3)),
    hjust = 1,
    vjust = -0.5
  ) +
  theme_bw()

dir.create("outputs", showWarnings = FALSE)
ggsave(
  "outputs/07a_edfig1_observed_vs_rarefied.pdf",
  p,
  width = 7,
  height = 6
)
ggsave(
  "outputs/07a_edfig1_observed_vs_rarefied.png",
  p,
  width = 7,
  height = 6,
  dpi = 300
)
