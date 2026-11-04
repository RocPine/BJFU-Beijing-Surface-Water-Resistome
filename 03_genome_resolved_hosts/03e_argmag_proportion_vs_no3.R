# Fig 3e: Proportion of ARG-carrying MAGs vs NO3-N concentration
# 数据: SI Table 8 (MAG quality) + Table 11 (MAG-level ARG carriage) + Table 7 (NO3_N)
# x = NO3-N 浓度, y = 携带 ARG 的 MAGs 占该样本所有 MAGs 的比例 (%)
# 运行: Rscript 03_genome_resolved_hosts/03e_argmag_proportion_vs_no3.R
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(dplyr)

norm_smp <- function(s) {
  sub("_binning.*$", "", as.character(s))
}

# Table 11: MAG 级 ARG 携带
t10 <- read.xlsx("data/Supplementary_Tables.xlsx", sheet = "Table 11",
                 startRow = 3, check.names = FALSE)
t10$site <- norm_smp(t10[, 1])
arg_mag_n <- t10 %>% distinct(site, genome = t10[, 1]) %>% count(site, name = "n_arg_mag")

# Table 8: 每样本 MAG 总数
t8 <- read.xlsx("data/Supplementary_Tables.xlsx", sheet = "Table 8",
                startRow = 2, check.names = FALSE)
t8$site <- norm_smp(t8$Sample)
mag_total <- t8 %>% count(site, name = "n_mag_total")

# Table 7: NO3-N
t7 <- read.xlsx("data/Supplementary_Tables.xlsx", sheet = "Table 7",
                startRow = 2, check.names = FALSE)
env <- t7[, c("Sample", "NO3_N")]
names(env)[1] <- "site"

data <- env %>%
  left_join(mag_total, by = "site") %>%
  left_join(arg_mag_n, by = "site") %>%
  mutate(n_arg_mag = ifelse(is.na(n_arg_mag), 0, n_arg_mag)) %>%
  filter(!is.na(NO3_N), !is.na(n_mag_total)) %>%
  mutate(arg_mag_pct = n_arg_mag / n_mag_total * 100)

r <- cor.test(data$NO3_N, data$arg_mag_pct, method = "spearman")
cat(sprintf("Fig 3e: rho=%.3f P=%.2e n=%d\n", r$estimate, r$p.value, nrow(data)))

p <- ggscatter(data, x = "NO3_N", y = "arg_mag_pct",
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
    y = "Proportion of ARG-carrying MAGs per sample (%)",
    x = "NO3-N concentration"
  )

if (!dir.exists("outputs")) dir.create("outputs")
ggsave("outputs/03e_hostabund_vs_no3.pdf", p, width = 6, height = 5)
ggsave("outputs/03e_hostabund_vs_no3.png", p, width = 6, height = 5, dpi = 300)
