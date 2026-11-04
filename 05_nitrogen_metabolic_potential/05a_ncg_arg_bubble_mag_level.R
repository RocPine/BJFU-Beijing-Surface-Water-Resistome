# Fig 5a (MAG 级): NCG-ARG 共宿主气泡图
# 数据: SI Table 11 (MAG 携带 ARG) + 05a_bubble_data.xlsx (df_count / N_carrying_MAGs)
# 注: df_ARG (MAG 携带 ARG) 来自 SI Table 11; MAG 级 N 基因携带 (df_N) 与
#     物种基因计数 (df_count) 无对应 SI 表, 读取 05a_bubble_data.xlsx
# 运行: Rscript 05_nitrogen_metabolic_potential/05a_ncg_arg_bubble_mag_level.R
library(tidyverse)
library(readxl)
library(ggh4x)
library(openxlsx)

df_ARG <- read.xlsx("data/Supplementary_Tables.xlsx", sheet = "Table 11",
                    startRow = 3, check.names = FALSE) %>%
  select(1) %>%
  distinct() %>%
  rename(genome = 1) %>%
  mutate(has_ARG = TRUE)

df_N <- read_excel(
  "05_nitrogen_metabolic_potential/05a_bubble_data.xlsx",
  sheet = "N_carrying_MAGs"
)

df_count <- read_excel(
  "05_nitrogen_metabolic_potential/05a_bubble_data.xlsx",
  sheet = "df_count"
)

top_species <- df_count %>%
  group_by(species) %>%
  summarise(total = sum(value), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 50) %>%
  pull(species)

species_mag_total <- df_N %>%
  filter(species %in% top_species) %>%
  distinct(species, genome) %>%
  count(species, name = "n_mag_total")

gene_mag_count <- df_N %>%
  filter(species %in% top_species) %>%
  distinct(species, gene, genome) %>%
  count(species, gene, name = "n_mag_gene")

df_bubble_size <- gene_mag_count %>%
  left_join(species_mag_total, by = "species") %>%
  mutate(prop_mag = n_mag_gene / n_mag_total)

df_N_ARG <- df_N %>%
  left_join(df_ARG, by = "genome") %>%
  mutate(has_ARG = ifelse(is.na(has_ARG), FALSE, TRUE))

df_star <- df_N_ARG %>%
  filter(species %in% top_species, has_ARG) %>%
  distinct(species, gene) %>%
  mutate(star = TRUE)

df_plot <- df_count %>%
  filter(species %in% top_species) %>%
  left_join(df_bubble_size, by = c("species", "gene")) %>%
  left_join(df_star, by = c("species", "gene")) %>%
  mutate(
    prop_mag = replace_na(prop_mag, 0),
    star = replace_na(star, FALSE)
  )

gene_pathway <- tribble(
  ~gene,  ~pathway,
  "amoa", "Nitrification",
  "hao",  "Nitrification",
  "hcp",  "Nitrification",
  
  "nirk", "Denitrification",
  "nirs", "Denitrification",
  "norb", "Denitrification",
  "norc", "Denitrification",
  "nosz", "Denitrification",
  "narg", "Denitrification",
  "narh", "Denitrification",
  "nari", "Denitrification",
  "narj", "Denitrification",
  "nary", "Denitrification",
  "narz", "Denitrification",
  "narc", "Denitrification",
  "napa", "Denitrification",
  "napb", "Denitrification",
  "napc", "Denitrification",
  
  "nifh", "N-fixation",
  "nifd", "N-fixation",
  "nifk", "N-fixation",
  "nifw", "N-fixation",
  
  "nasa", "Assimilatory nitrate reduction",
  "nasb", "Assimilatory nitrate reduction",
  "nira", "Assimilatory nitrate reduction",
  "narb", "Assimilatory nitrate reduction",
  
  "nirb", "Dissimilatory nitrate reduction",
  "nird", "Dissimilatory nitrate reduction",
  "nrfa", "Dissimilatory nitrate reduction",
  "nrfb", "Dissimilatory nitrate reduction",
  "nrfc", "Dissimilatory nitrate reduction",
  
  "urea", "Other nitrogen metabolism",
  "ureb", "Other nitrogen metabolism",
  "urec", "Other nitrogen metabolism",
  "glna", "Other nitrogen metabolism",
  "gs",   "Other nitrogen metabolism",
  "gdh",  "Other nitrogen metabolism",
  "ansb", "Other nitrogen metabolism",
  "nmo",  "Other nitrogen metabolism"
)

pathway_order <- c(
  "Nitrification",
  "Denitrification",
  "N-fixation",
  "Assimilatory nitrate reduction",
  "Dissimilatory nitrate reduction",
  "Other nitrogen metabolism"
)

df_plot <- df_plot %>%
  left_join(gene_pathway, by = "gene") %>%
  mutate(
    pathway = factor(pathway, levels = pathway_order)
  ) %>%
  arrange(pathway, gene) %>%
  mutate(
    gene = factor(gene, levels = unique(gene))
  )

x_levels <- c("Pathway", unique(as.character(df_plot$species)))

df_plot <- df_plot %>%
  mutate(
    species_num = as.numeric(factor(species, levels = unique(species))) + 1
  )

gene_order <- df_plot %>%
  arrange(pathway, gene) %>%
  pull(gene) %>%
  unique()

df_plot <- df_plot %>%
  mutate(gene = factor(gene, levels = gene_order))

p <- ggplot() +
  geom_tile(
    data = df_plot %>% distinct(gene, pathway),
    aes(x = 1, y = gene, fill = pathway),
    width = 0.8,
    height = 1.0
  ) +
  geom_point(
    data = df_plot,
    aes(x = species_num, y = gene, size = prop_mag, color = value),
    alpha = 0.85
  ) +
  geom_text(
    data = df_plot %>% filter(star),
    aes(x = species_num, y = gene, label = "*"),
    color = "#FF2D2D",
    size = 5.5,
    fontface = "bold"
  ) +
  scale_x_continuous(
    breaks = 1:length(x_levels),
    labels = x_levels,
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_size_continuous(
    name = "Proportion of MAGs\ncarrying this gene",
    range = c(1.2, 6)
  ) +
  scale_color_gradient(
    low = "#C6DBEF",
    high = "#08306B",
    trans = "sqrt",
    name = "Gene count"
  ) +
  scale_fill_manual(
    name = "Nitrogen cycling pathway",
    values = c(
      "Nitrification" = "#1F78B4",
      "Denitrification" = "#33A02C",
      "N-fixation" = "#6A3D9A",
      "Assimilatory nitrate reduction" = "#FF7F00",
      "Dissimilatory nitrate reduction" = "#E31A1C",
      "Other nitrogen metabolism" = "#BDBDBD"
    )
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
    axis.text.y = element_text(size = 8),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )

p

if (!dir.exists("outputs")) dir.create("outputs")
ggsave("outputs/05a_ncg_arg_bubble_mag_level.pdf", p, width = 14, height = 8)
ggsave("outputs/05a_ncg_arg_bubble_mag_level.png", p, width = 14, height = 8, dpi = 300)
