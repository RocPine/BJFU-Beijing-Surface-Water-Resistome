# Fig 1b: ARG class composition circle packing plot
# 数据: SI Table 4 (Sample x ARG subtype; 列名格式 "class__subtype")
# 运行: Rscript 01_spatial_patterns/01b_arg_class_circle_packing.R
library(openxlsx)
library(dplyr)
library(tidyr)
library(tidygraph)
library(ggraph)

df <- read.xlsx("data/Supplementary_Tables.xlsx", sheet = "Table 4",
                startRow = 2, check.names = FALSE)

subtype_cols <- names(df)[-1]
class_of <- sub("__.*$", "", subtype_cols)
subtype_of <- sub("^[^_]+__", "", subtype_cols)
value <- colSums(df[, -1], na.rm = TRUE)

df <- data.frame(Type = class_of, Subtype = subtype_of, Value = value,
                 stringsAsFactors = FALSE) %>%
  filter(Value > 0)

nodes_type <- df %>%
  group_by(name = Type) %>%
  summarise(size = sum(Value, na.rm = TRUE), .groups = "drop") %>%
  mutate(Type = name)
nodes_sub <- df %>%
  rename(name = Subtype) %>%
  select(name, size = Value, Type)

vertices <- bind_rows(nodes_type, nodes_sub)
edges <- df %>% select(from = Type, to = Subtype)
graph <- tbl_graph(nodes = vertices, edges = edges, directed = TRUE)

top_sub <- df %>%
  group_by(Type) %>%
  group_modify(~{
    n_sub <- nrow(.x)
    n_top <- ifelse(n_sub <= 5, 2, 10)
    slice_max(.x, order_by = Value, n = n_top)
  }) %>%
  pull(Subtype)

macaron_30 <- c(
  "#FFC8A2", "#FFB3BA", "#E0BBE4", "#FFDFBA", "#FFF9C4",
  "#FFFFBA", "#F0FFF0", "#BAFFC9", "#E0F2F1", "#BAE1FF",
  "#E3F2FD", "#FEE5D9", "#F3E5F5", "#FADADD", "#FFF0F5",
  "#F9EB97", "#F5F5DC", "#D4F0F0", "#FAF0E6", "#D9E4FF",
  "#85C1E9", "#FF96AD", "#E6E6FA", "#F0FFFF", "#CCCCFF",
  "#E2F0CB", "#F5FFFA", "#B5EAD7", "#FFFACD", "#C7CEEA",
  "#FF9AA2", "#FFB347", "#FFFF99", "#B5EAEA", "#A0CED9",
  "#C2F0F0", "#D5AAFF", "#FFC0CB"
)

n_type <- length(unique(df$Type))
if (n_type > length(macaron_30)) {
  macaron_30 <- c(macaron_30, grDevices::hcl.colors(n_type - length(macaron_30), "Pastel 1"))
}

p <- ggraph(graph, layout = 'circlepack', weight = size) +
  geom_node_circle(aes(filter = leaf, fill = Type), color = "#E0E0E0", linewidth = 0.2) +
  geom_node_circle(aes(filter = !leaf), fill = NA, color = "grey", linewidth = 0.5) +
  geom_node_text(aes(label = name, filter = name %in% top_sub),
                 size = 3, check_overlap = TRUE) +
  scale_fill_manual(values = macaron_30, name = "Type") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  coord_fixed()

if (!dir.exists("outputs")) dir.create("outputs")
ggsave("outputs/01b_arg_class_circle_packing.pdf", p, width = 8, height = 6)
ggsave("outputs/01b_arg_class_circle_packing.png", p, width = 8, height = 6, dpi = 300)
