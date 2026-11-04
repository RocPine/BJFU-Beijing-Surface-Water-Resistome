# Fig 5c: NCG abundance differences between ARG hosts and non-hosts (volcano)
# 数据: SI Table 16 (45-NCG Wilcoxon: Gene, mean_Host, mean_Nonhost, logFC, pval, padj, Function, group)
# 运行: Rscript 05_nitrogen_metabolic_potential/05c_ncg_volcano.R
library(openxlsx)
library(ggplot2)
library(ggrepel)

res <- read.xlsx("data/Supplementary_Tables.xlsx", sheet = "Table 16",
                 startRow = 2, check.names = FALSE)

if (!"group" %in% names(res)) {
  res$group <- "NotSig"
  res$group[res$padj < 0.05 & res$logFC > 0.585] <- "Up"
  res$group[res$padj < 0.05 & res$logFC < -0.585] <- "Down"
}

sig_df <- subset(res, group != "NotSig")

p <- ggplot(res, aes(x = logFC, y = -log10(padj))) +
  geom_point(
    data = subset(res, group == "NotSig"),
    color = "grey", size = 2, alpha = 0.6
  ) +
  geom_point(
    data = sig_df,
    aes(color = Function), size = 3.5, alpha = 0.9
  ) +
  scale_color_discrete() +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 13),
    axis.line = element_line(color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  ) +
  labs(
    x = "log2 Fold Change (Host vs Non-host)",
    y = "-log10(FDR)",
    color = "Function"
  ) +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = sig_df,
    aes(label = Gene, color = Function),
    size = 4, max.overlaps = 50, show.legend = FALSE
  )

if (!dir.exists("outputs")) dir.create("outputs")
ggsave("outputs/05c_ncg_volcano.pdf", p, width = 8, height = 6)
ggsave("outputs/05c_ncg_volcano.png", p, width = 8, height = 6, dpi = 300)
