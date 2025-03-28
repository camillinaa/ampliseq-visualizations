# Import libraries

library(data.tree)
library(ggtree)
library(ggplot2)
library(ape)
library(dplyr)
library(stringr)


dir = "rectal-cancer"

lfc <- read_delim(paste0("input/",dir,"/ancombc-formula/lfc_slice-level6.csv"), delim=",")
qval <- read_delim(paste0("input/",dir,"/ancombc-formula/q_val_slice-level6.csv"), delim=",") 

# Parse taxonomy

taxonomy_clean <- lfc %>%
  select(id) %>%
  distinct() %>%
  filter(str_detect(id, ";")) %>%  # Remove malformed or missing taxonomies
  separate(id, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";", fill = "right") %>%
  mutate(across(everything(), ~ ifelse(.x == "", NA, .x)))

# Create edges

levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

edges <- map2_df(levels[-length(levels)], levels[-1], ~ {
  taxonomy_clean %>%
    filter(!is.na(.data[[.x]]), !is.na(.data[[.y]])) %>%
    distinct(parent = .data[[.x]], child = .data[[.y]]) %>%
    mutate(level = .y)
})

# Build tree using ape

edges$pathString <- paste("Life", edges$parent, edges$child, sep = "/") ### NOT OK
taxonomy_tree <- as.Node(edges)
phylo_tree <- as.phylo(taxonomy_tree)

# Annotate tree

data_agg <- merge(lfc, qval, by = "id", suffixes = c("_lfc", "_qval"))

data_agg <- data_agg %>%
  separate(id, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";", fill = "right") %>%
  unite("tax_label", Kingdom:Genus, sep = ";", remove = FALSE, na.rm = TRUE) %>%
  mutate(
    most_specific = coalesce(Genus, Family, Order, Class, Phylum, Kingdom),  # fallback to most specific
    point_size = abs(conditioncancer_lfc),
    lfc_color = case_when(
      conditioncancer_lfc > 0 ~ "Up",
      conditioncancer_lfc < 0 ~ "Down",
      TRUE ~ "Neutral"
    ),
    lfc_color = factor(lfc_color, levels = c("Up", "Down", "Neutral")),
    sig_shape = ifelse(conditioncancer_qval < 0.05, 17, 16)
  )


# Create tree

p <- ggtree(phylo_tree, layout = "circular", branch.length = "none")

p$data <- p$data %>%
  left_join(data_agg, by = c("label" = "most_specific"))  # Join on most specific taxon

p +
  geom_point(aes(color = lfc_color, size = point_size, shape = factor(sig_shape)), na.rm = TRUE) +
  geom_tiplab(aes(label = label), size = 2.8, hjust = -0.15) +
  geom_label2(
    data = subset(p$data, isTip == FALSE),
    aes(label = label), 
    hjust = -0.4,
    size = 2.8,
    color = "darkgreen"
  ) +
  scale_color_manual(
    values = c("Up" = "red", "Down" = "blue", "Neutral" = "grey"),
    name = "LFC Direction"
  ) +
  scale_size_continuous(
    range = c(2, 6),
    name = expression("|log"[2]*" fold change|")
  ) +
  scale_shape_manual(
    values = c("16" = 16, "17" = 17),
    labels = c("Not significant", "q < 0.05"),
    name = "Significance"
  ) +
  labs(title = "Genus-Level Cladogram of Cancer-Associated Microbiome Changes (ANCOM-BC)") +
  theme(legend.position = "right")


ggsave("results/ancombc/cladogram-level6.png", plot = last_plot())
