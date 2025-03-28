library(data.tree)
library(ggtree)
library(ggplot2)
library(ape)
library(dplyr)
library(stringr)
library(readr)
library(purrr)

dir <- "rectal-cancer"


# --- BUILD TREE STRUCTURE --- #


# Use most specific taxonomic rank (6 - genus)

lfc <- read_delim(paste0("input/",dir,"/ancombc-formula/lfc_slice-level6.csv"), delim=",")
qval <- read_delim(paste0("input/",dir,"/ancombc-formula/q_val_slice-level6.csv"), delim=",") 

# Parse taxonomy

taxonomy <- lfc %>%
  select(id) %>%
  distinct() %>%
  filter(str_detect(id, ";")) %>%  # Remove malformed or missing taxonomies
  separate(id, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";", fill = "right") %>%
  mutate(across(everything(), ~ ifelse(.x == "", NA, .x)))

levels <- colnames(taxonomy)

taxonomy_suff <- taxonomy %>%
  mutate(across(all_of(levels), 
                .fns = ~ ifelse(is.na(.x), NA, paste0(.x, "__", tolower(substr(cur_column(), 1, 1))))))

taxonomy_suff$pathString <- apply(taxonomy_suff, 1, function(row) {
  paste(c("Life", row[!is.na(row)]), collapse = "/")
}) # add pathstring for as.Node()

taxonomy_tree <- as.Node(taxonomy_suff)
phylo_tree <- as.phylo(taxonomy_tree)


# --- BUILD ANNOTATION USING ALL FILES --- #


level_range <- 2:6  # Levels from Phylum (2) to Genus (6)

load_level_data <- function(level) {
  lfc_path <- paste0("input/", dir, "/ancombc-formula/lfc_slice-level", level, ".csv")
  qval_path <- paste0("input/", dir, "/ancombc-formula/q_val_slice-level", level, ".csv")
  
  lfc_df <- read_delim(lfc_path, delim = ",")
  qval_df <- read_delim(qval_path, delim = ",")
  
  df <- merge(lfc_df, qval_df, by = "id", suffixes = c("_lfc", "_qval"))
  
  tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
  
  df <- df %>%
    filter(str_detect(id, ";")) %>%
    separate(id, into = tax_levels, sep = ";", fill = "right", remove = FALSE) %>%
    mutate(
      tax_level = tax_levels[level],
      label = paste(.data[[tax_levels[level]]], tolower(substring(tax_level, 1, 1)), sep = "__")
    ) %>%
    filter(!is.na(label)) %>%
    mutate(
      point_size = abs(conditioncancer_lfc),
      lfc_color = case_when(
        conditioncancer_lfc > 0 ~ "Up",
        conditioncancer_lfc < 0 ~ "Down",
        TRUE ~ "Neutral"
      ),
      lfc_color = factor(lfc_color, levels = c("Up", "Down", "Neutral")),
      sig_shape = ifelse(conditioncancer_qval < 0.05, 17, 16)
    ) %>%
    select(label, tax_level, point_size, lfc_color, sig_shape)
  
  return(df)
}


# Combine data from all levels
data_agg_all_levels <- map_dfr(level_range, load_level_data)


# --- CREATE CLADOGRAM ---#


p <- ggtree(phylo_tree, layout = "circular", branch.length = "none")

p$data <- p$data %>%
  left_join(data_agg_all_levels, by = "label") %>%
  mutate(
    label_clean = gsub("__.*", "", label)  # Create a new column with the cleaned label
  )


p +
  geom_point(aes(color = lfc_color, size = point_size, shape = factor(sig_shape)), na.rm = TRUE) +
  geom_tiplab(aes(label = label_clean), size = 2.8, hjust = -0.15) +
  geom_label2(
    data = subset(p$data, isTip == FALSE),
    aes(label = label_clean), 
    hjust = 0,
    size = 1.5,
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
  labs(title = "Genus-Level Cladogram of Condition/Treatment-Associated Microbiome Changes (ANCOM-BC)") +
  theme(legend.position = "right")

ggsave("results/ancombc/cladogram-main.png", plot = last_plot(), dpi = 600)


