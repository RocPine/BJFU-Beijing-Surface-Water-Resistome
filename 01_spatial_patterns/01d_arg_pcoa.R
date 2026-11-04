# Fig 1d: PCoA of ARG profiles (Bray-Curtis)
# 数据: SI Table 4 (Sample x ARG subtype abundance matrix)
# 运行: Rscript 01_spatial_patterns/01d_arg_pcoa.R
library(vegan)
library(ggplot2)
library(openxlsx)

otu <- read.xlsx("data/Supplementary_Tables.xlsx", sheet = "Table 4",
                 startRow = 2, check.names = FALSE)
rownames(otu) <- otu[, 1]
otu <- otu[, -1]
otu <- as.matrix(otu)
storage.mode(otu) <- "numeric"

otu.dist <- vegdist(otu, method = "bray")
pcoa <- cmdscale(otu.dist, eig = TRUE, k = 2)
pc12 <- as.data.frame(pcoa$points)
pc12$samples <- rownames(pc12)

eig <- pcoa$eig
pc_importance <- round(eig / sum(eig[eig > 0]) * 100, 2)

p <- ggplot(pc12, aes(x = V1, y = V2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey80") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
  geom_hex(bins = 30) +
  scale_fill_viridis_c(option = "C", name = "Density") +
  geom_density_2d(color = "#BDBDBD", size = 0.5, alpha = 0.7) +
  geom_point(size = 1, alpha = 0.7, color = "#66BB6A") +
  labs(
    x = paste0("PCoA1 (", pc_importance[1], "%)"),
    y = paste0("PCoA2 (", pc_importance[2], "%)")
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    axis.line = element_blank(),
    plot.title = element_text(size = 18, hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

if (!dir.exists("outputs")) dir.create("outputs")
ggsave("outputs/01d_arg_pcoa.pdf", p, width = 8, height = 6)
ggsave("outputs/01d_arg_pcoa.png", p, width = 8, height = 6, dpi = 300)
