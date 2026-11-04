# Fig 2a: RDA of ARG profiles against environmental variables
# 数据: SI Table 4 (ARG Type/Class abundance) + Table 7 (environmental variables)
# Fig. 2a, left: ARG Type/Class RDA.
# Run from repository root: Rscript 02_environmental_drivers/fig2a_arg_rda.R
library(vegan)
library(ggplot2)
library(openxlsx)
library(grid)

scale_factor <- (90 / 25.4) / 10

otu_matrix <- read.xlsx("data/Supplementary_Tables.xlsx", sheet = "Table 4",
                        startRow = 2, check.names = FALSE)
rownames(otu_matrix) <- otu_matrix[, 1]
otu_matrix <- otu_matrix[, -1]
otu_matrix <- as.matrix(otu_matrix)
storage.mode(otu_matrix) <- "numeric"
# Aggregate ARG subtype columns to ARG Type/Class using the text before "__".
arg_types <- sub("__.*$", "", colnames(otu_matrix))
otu_matrix <- sapply(unique(arg_types), function(x) {
  rowSums(otu_matrix[, arg_types == x, drop = FALSE])
})
colnames(otu_matrix) <- unique(arg_types)
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
rda1_pct <- 100 * xxxx$concont$importance["Proportion Explained", "RDA1"]
rda2_pct <- 100 * xxxx$concont$importance["Proportion Explained", "RDA2"]
constrained12_pct <- rda1_pct + rda2_pct
cat(sprintf("Constrained RDA axes: RDA1 = %.2f%%, RDA2 = %.2f%%, first two = %.2f%%\n",
            rda1_pct, rda2_pct, constrained12_pct))

pdat <- res$CCA
samples <- data.frame(sample = row.names(pdat$u), RDA1 = pdat$u[, 1], RDA2 = pdat$u[, 2])
species <- data.frame(spece = row.names(pdat$v), RDA1 = pdat$v[, 1], RDA2 = pdat$v[, 2])
envi <- data.frame(en = row.names(pdat$biplot), RDA1 = pdat$biplot[, 1], RDA2 = pdat$biplot[, 2])

p <- ggplot() +
  geom_point(data = samples, aes(x = RDA1, y = RDA2), size = 7 * scale_factor * 0.8, shape = 21, fill = '#B5ABC760', colour = '#545454', stroke = 0.2 * scale_factor * 0.8) +
  geom_segment(data = envi, aes(x = 0, xend = RDA1, y = 0, yend = RDA2), arrow = arrow(length = unit(0.5 * scale_factor, "cm")), colour = '#86B7DB', linewidth = 1.2 * scale_factor) +
  geom_text(data = envi, aes(x = RDA1 + 0.02, y = RDA2 + 0.05, label = en), size = 5.6 / 2.845276, colour = 'black', family = "Arial", check_overlap = FALSE) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(sprintf('RDA1 (%.2f%%)', rda1_pct)) + ylab(sprintf('RDA2 (%.2f%%)', rda2_pct)) +
  labs(title = "") +
  theme(axis.title.x = element_text(size = 7, family = "Arial", color = "black"),
        axis.title.y = element_text(size = 7, family = "Arial", angle = 90, color = "black"),
        axis.text.y = element_text(size = 7, family = "Arial", color = "black"),
        axis.text.x = element_text(size = 7, family = "Arial", color = "black"),
        panel.grid = element_blank(),
        text = element_text(family = "Arial"),
        legend.text = element_text(size = 7, family = "Arial"),
        plot.title = element_text(size = 7, family = "Arial"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1 * scale_factor)) +
  guides(color = guide_legend(title = NULL))

if (!dir.exists("outputs")) dir.create("outputs")
ggsave("outputs/02a_arg_rda.pdf", p, device = cairo_pdf, width = 90 / 25.4, height = 72 / 25.4, units = "in")
ggsave("outputs/02a_arg_rda.png", p, width = 90 / 25.4, height = 72 / 25.4, units = "in", dpi = 300)
