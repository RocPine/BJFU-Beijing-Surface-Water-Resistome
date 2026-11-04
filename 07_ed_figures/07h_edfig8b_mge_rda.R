# ED Fig 8b: Contributions of environmental variables to variation in MGE profiles (RDA)
# 数据: SI Table 14 (Sample x MGE subtype abundance) + Table 7 (environmental variables)
# 运行: Rscript 07_ed_figures/07h_edfig8b_mge_rda.R
library(vegan)
library(ggplot2)
library(openxlsx)

otu_matrix <- read.xlsx("data/Supplementary_Tables.xlsx", sheet = "Table 14",
                        startRow = 2, check.names = FALSE)
rownames(otu_matrix) <- otu_matrix[, 1]
otu_matrix <- otu_matrix[, -1]
otu_matrix <- as.matrix(otu_matrix)
storage.mode(otu_matrix) <- "numeric"
otu_matrix <- otu_matrix[rowSums(otu_matrix) > 0, colSums(otu_matrix) > 0]

envdat_raw <- read.xlsx("data/Supplementary_Tables.xlsx", sheet = "Table 7",
                        startRow = 2, check.names = FALSE)
rownames(envdat_raw) <- envdat_raw[, 1]
envdat_raw <- envdat_raw[, -1]
envdat_raw <- envdat_raw[, sapply(envdat_raw, is.numeric)]
envdat <- envdat_raw[match(rownames(otu_matrix), rownames(envdat_raw)), ]
envdat <- envdat[, colSums(is.na(envdat)) == 0]
keep <- complete.cases(envdat)
otu_matrix <- otu_matrix[keep, ]
envdat <- envdat[keep, ]

decorana(otu_matrix)
res <- rda(otu_matrix ~ ., envdat)
xxxx <- summary(res)

aa <- xxxx$concont$importance
aa <- round(aa, 4)

pdat <- res$CCA
samples <- data.frame(sample = row.names(pdat$u), RDA1 = pdat$u[, 1], RDA2 = pdat$u[, 2])
envi <- data.frame(en = row.names(pdat$biplot), RDA1 = pdat$biplot[, 1], RDA2 = pdat$biplot[, 2])

p <- ggplot() +
  geom_hline(aes(yintercept = 0), colour = "gray88", linetype = "dashed") +
  geom_vline(aes(xintercept = 0), colour = "gray88", linetype = "dashed") +
  geom_point(data = samples, aes(x = RDA1, y = RDA2), size = 7, shape = 21, fill = '#B5ABC760', colour = '#545454', stroke = 0.2) +
  geom_segment(data = envi, aes(x = 0, xend = RDA1, y = 0, yend = RDA2), arrow = arrow(length = unit(0.5, "cm")), colour = '#86B7DB', size = 1.2) +
  geom_text(data = envi, aes(x = RDA1 + 0.02, y = RDA2 + 0.05, label = en), size = 10, colour = 'black', family = "serif", check_overlap = TRUE) +
  theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) +
  xlab(paste('RDA1 (', aa[2, 1] * 100, '%)', sep = '')) + ylab(paste('RDA2 (', aa[2, 2] * 100, '%)', sep = '')) +
  labs(title = "") +
  theme(axis.title.x = element_text(size = 35, family = "serif", color = "black"),
        axis.title.y = element_text(size = 35, family = "serif", angle = 90, color = "black"),
        axis.text.y = element_text(size = 30, family = "serif", color = "black"),
        axis.text.x = element_text(size = 30, family = "serif", color = "black"),
        panel.grid = element_blank(),
        text = element_text(family = "serif"),
        legend.text = element_text(size = 20, family = "serif"),
        plot.title = element_text(size = 35, family = "serif"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  guides(color = guide_legend(title = NULL))

if (!dir.exists("outputs")) dir.create("outputs")
ggsave("outputs/07h_edfig8b_mge_rda.pdf", p, width = 10, height = 8)
ggsave("outputs/07h_edfig8b_mge_rda.png", p, width = 10, height = 8, dpi = 300)
