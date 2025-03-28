library(data.tree)
library(ggtree)
library(ggplot2)
library(ape)
library(dplyr)
library(stringr)

dir = "rectal-cancer"

lfc <- read_delim(paste0("input/",dir,"/ancombc-formula/lfc_slice-2.csv"), delim=",")
qval <- read_delim(paste0("input/",dir,"/ancombc-formula/q_val_slice-2.csv"), delim=",")

# Parse taxonomy

taxonomy_clean <- lfc %>%
  select(id) %>%
  distinct() %>%
  filter(str_detect(id, ";")) %>%  # Remove malformed or missing taxonomies
  separate(id, into = c("Kingdom", "Phylum"), sep = ";", fill = "right") %>%
  filter(!is.na(Kingdom) & !is.na(Phylum)) %>%
  mutate(Kingdom = ifelse(Kingdom == "", "Bacteria", Kingdom))  # Default Kingdom

# Create edges

edges <- taxonomy_clean %>%
  distinct(Kingdom, Phylum) %>%
  rename(parent = Kingdom, child = Phylum)

# Build tree using ape

edges$pathString <- paste("Life", edges$parent, edges$child, sep = "/")
taxonomy_tree <- as.Node(edges)
phylo_tree <- as.phylo(taxonomy_tree) # convert to phylo object

# Annotate tree

data_agg <- merge(lfc, qval, by = "id", suffixes = c("_lfc", "_qval"))

data_agg <- data_agg %>%
  mutate(
    # Extract Phylum-level label
    Phylum = str_split_fixed(id, ";", 2)[, 2],
    
    # Color by LFC direction
    lfc_color = case_when(
      conditioncancer_lfc > 0 ~ "Up",
      conditioncancer_lfc < 0 ~ "Down",
      TRUE ~ "Neutral"
    ),
    lfc_color = factor(lfc_color, levels = c("Up", "Down", "Neutral")),
    
    # Size: magnitude of log-fold change
    point_size = abs(conditioncancer_lfc),
    
    # Shape: use star (8) if significant, else circle (16)
    sig_shape = ifelse(conditioncancer_qval < 0.05, 8, 16)
  )


p <- ggtree(phylo_tree, layout = "circular", branch.length = "none")

p$data <- p$data %>%
  left_join(data_agg, by = c("label" = "Phylum"))  # label is the tip name

p +
  geom_point(aes(color = lfc_color, size = point_size, shape = factor(sig_shape)), na.rm = TRUE) +
  geom_tiplab(aes(label = label), size = 3) +
  scale_color_manual(
    values = c("Up" = "red", "Down" = "blue", "Neutral" = "grey"),
    name = "LFC Direction"
  ) +
  scale_size_continuous(
    range = c(2, 6),
    name = expression("|log"[2]*" fold change|")
  ) +
  scale_shape_manual(
    values = c("16" = 16, "8" = 8),
    labels = c("Not significant", "q < 0.05"),
    name = "Significance"
  ) +
  labs(
    title = "Differentially Abundant Taxa in Condition/Treatment via ANCOM-BC"
  ) +
  theme(legend.position = "right")

