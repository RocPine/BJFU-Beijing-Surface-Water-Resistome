# Fig 3d: Number of MAG-level ARG subtypes vs NO3-N concentration
# 数据: SI Table 8 (MAG quality) + Table 11 (MAG-level ARG carriage) + Table 7 (NO3_N)
# x = NO3-N 浓度, y = 每样本 MAG 水平检出的去重 ARG 亚型数
# 注: 与 Fig 3e (ARG-carrying MAG 占比 vs NO3-N) 互补: 3d 反映耐药基因丰富度, 3e 反映携带谱系占比
# 运行: Rscript 03_genome_resolved_hosts/03d_argsubtypes_vs_no3.R
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(dplyr)

# 归一化样本名 (去 _binning 后缀)
norm_smp <- function(s) {
  sub("_binning.*$", "", as.character(s))
}

# Table 11: MAG 级 ARG 携带 (genome, ARG)
t10 <- read.xlsx("data/Supplementary_Tables.xlsx", sheet = "Table 11",
                 startRow = 3, check.names = FALSE)
t10$site <- norm_smp(t10[, 1])
t10$arg <- as.character(t10[, 2])

# 每样本: 去重 ARG 亚型数
arg_sub_n <- t10 %>% distinct(site, arg) %>% count(site, name = "n_arg_sub")

# Table 7: NO3-N
t7 <- read.xlsx("data/Supplementary_Tables.xlsx", sheet = "Table 7",
                startRow = 2, check.names = FALSE)
env <- t7[, c("Sample", "NO3_N")]
names(env)[1] <- "site"

data <- env %>%
  left_join(arg_sub_n, by = "site") %>%
  mutate(n_arg_sub = ifelse(is.na(n_arg_sub), 0, n_arg_sub)) %>%
  filter(!is.na(NO3_N))

r <- cor.test(data$NO3_N, data$n_arg_sub, method = "spearman")
cat(sprintf("Fig 3d: rho=%.3f P=%.2e n=%d\n", r$estimate, r$p.value, nrow(data)))

p <- ggscatter(data, x = "NO3_N", y = "n_arg_sub",
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
    y = "Number of ARG subtypes carried by MAGs per sample",
    x = "NO3-N concentration"
  )

if (!dir.exists("outputs")) dir.create("outputs")
ggsave("outputs/03d_argnum_vs_hostabund.pdf", p, width = 6, height = 5)
ggsave("outputs/03d_argnum_vs_hostabund.png", p, width = 6, height = 5, dpi = 300)
