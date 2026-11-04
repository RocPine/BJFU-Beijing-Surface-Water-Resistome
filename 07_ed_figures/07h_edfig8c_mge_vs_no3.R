# ED Fig 8c: Total MGE abundance vs NO3-N concentration (scatter + Spearman)
# 数据: SI Table 14 (MGE abundance) + Table 7 (NO3_N)
# 运行: Rscript 07_ed_figures/07h_edfig8c_mge_vs_no3.R
library(openxlsx)
library(ggplot2)
library(ggpubr)

t14 <- read.xlsx("data/Supplementary_Tables.xlsx", sheet = "Table 14",
                 startRow = 2, check.names = FALSE)
t7 <- read.xlsx("data/Supplementary_Tables.xlsx", sheet = "Table 7",
                startRow = 2, check.names = FALSE)

mge_mat <- as.matrix(t14[, -1])
storage.mode(mge_mat) <- "numeric"
mge_total <- data.frame(sample = t14[, 1], mge_abund = rowSums(mge_mat, na.rm = TRUE))
env <- t7[, c("Sample", "NO3_N")]

data <- merge(mge_total, env, by.x = "sample", by.y = "Sample")

p <- ggscatter(data, x = "NO3_N", y = "mge_abund",
               color = "#545454", fill = "#93B8DB", size = 2.5, shape = 21,
               stroke = 0.6, conf.int = TRUE,
               add.params = list(fill = "#93B8DB")) +
  stat_cor(method = "spearman", size = 4, fontface = "plain") +
  stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                        size = 4, fontface = "plain") +
  geom_smooth(method = "lm", color = "#65677F", linetype = "solid") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line.y.left = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16)
  ) +
  labs(
    y = "Total MGE abundance",
    x = "NO3--N concentration"
  )

if (!dir.exists("outputs")) dir.create("outputs")
ggsave("outputs/07h_edfig8c_mge_vs_no3.pdf", p, width = 6, height = 5)
ggsave("outputs/07h_edfig8c_mge_vs_no3.png", p, width = 6, height = 5, dpi = 300)
