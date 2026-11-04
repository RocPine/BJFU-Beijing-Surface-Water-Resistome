# Fig 1c: Phylum-level taxonomic composition circle packing plot
# 数据: SI Table 5 (Species x samples + 分类列; 门/属组成由 Table 5 聚合)
# 运行: Rscript 01_spatial_patterns/01c_phylum_circle_packing.R
library(openxlsx)
library(dplyr)
library(tidygraph)
library(ggraph)

df <- read.xlsx("data/Supplementary_Tables.xlsx", sheet = "Table 5",
                startRow = 2, check.names = FALSE)

abund <- df[, -c(1:7)]
is_bact <- df$Kingdom == "Bacteria"
bact_total <- colSums(abund[is_bact, , drop = FALSE], na.rm = TRUE)
pct_per_sample <- sweep(abund[is_bact, , drop = FALSE], 2, bact_total, "/") * 100
value <- rowMeans(pct_per_sample, na.rm = TRUE)

df <- data.frame(
  Phylum = df$Phylum[is_bact],
  Genus = df$Genus[is_bact],
  Value = value,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(Phylum), Phylum != "", Value > 0)

nodes_phylum <- df %>%
  group_by(name = Phylum) %>%
  summarise(size = sum(Value, na.rm = TRUE), .groups = "drop") %>%
  mutate(Phylum = name)
nodes_genus <- df %>%
  rename(name = Genus) %>%
  select(name, size = Value, Phylum)

vertices <- bind_rows(nodes_phylum, nodes_genus)
edges <- df %>% select(from = Phylum, to = Genus)
graph <- tbl_graph(nodes = vertices, edges = edges, directed = TRUE)

top_genus <- df %>%
  group_by(Phylum) %>%
  group_modify(~{
    n_sub <- nrow(.x)
    n_top <- ifelse(n_sub <= 5, 2, 10)
    slice_max(.x, order_by = Value, n = n_top)
  }) %>%
  pull(Genus)

macaron_43 <- c(
  "#FFC8A2", "#FFB3BA", "#E0BBE4", "#FFDFBA", "#FFF9C4",
  "#FFFFBA", "#F0FFF0", "#BAFFC9", "#E0F2F1", "#BAE1FF",
  "#E3F2FD", "#FEE5D9", "#F3E5F5", "#FADADD", "#FFF0F5",
  "#F9EB97", "#F5F5DC", "#D4F0F0", "#FAF0E6", "#D9E4FF",
  "#85C1E9", "#FF96AD", "#E6E6FA", "#F0FFFF", "#CCCCFF",
  "#E2F0CB", "#F5FFFA", "#B5EAD7", "#FFFACD", "#C7CEEA",
  "#FF9AA2", "#FFB347", "#FFFF99", "#B5EAEA", "#A0CED9",
  "#C2F0F0", "#D5AAFF", "#FFC0CB", "#FFE4B5", "#C1FFC1",
  "#ADD8E6", "#F5DEB3", "#E6E6FA"
)

n_phylum <- length(unique(df$Phylum))
if (n_phylum > length(macaron_43)) {
  macaron_43 <- c(macaron_43, grDevices::hcl.colors(n_phylum - length(macaron_43), "Pastel 1"))
}

p <- ggraph(graph, layout = 'circlepack', weight = size) +
  geom_node_circle(aes(filter = leaf, fill = Phylum), color = "#E0E0E0", linewidth = 0.2) +
  geom_node_circle(aes(filter = !leaf), fill = NA, color = "grey", linewidth = 0.5) +
  geom_node_text(aes(label = name, filter = name %in% top_genus),
                 size = 3, check_overlap = TRUE) +
  scale_fill_manual(values = macaron_43, name = "Phylum") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  coord_fixed()

if (!dir.exists("outputs")) dir.create("outputs")
ggsave("outputs/01c_phylum_circle_packing.pdf", p, width = 8, height = 6)
ggsave("outputs/01c_phylum_circle_packing.png", p, width = 8, height = 6, dpi = 300)
