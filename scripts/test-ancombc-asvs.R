library(ggtree)
library(ape)
library(dplyr)
library(ggplot2)

dir = "rectal-cancer"

# Step 1: Read the tree from Newick format
tree <- read.tree(paste0("input/",dir,"/ancombc-formula/tree.nwk")) # /results/qiime2/phylogenetic_tree/tree.nwk

# Step 2: Prepare your taxonomy table that links ASVs to taxonomic groups
asv_taxonomy <- read_tsv(paste0("input/",dir,"/ancombc-formula/ASV_tax.gtdb.tsv"))
asv_taxonomy[is.na(asv_taxonomy)] <- ""
asv_taxonomy$label <- paste0(asv_taxonomy$Kingdom, ";", asv_taxonomy$Phylum, ";", asv_taxonomy$Class, ";", asv_taxonomy$Order, ";", asv_taxonomy$Family, ";", asv_taxonomy$Genus)
asv_taxonomy <- asv_taxonomy[c("ASV_ID","label")]

# Step 3: Prepare your aggregated LFC and q-value data at taxonomic level
lfc <- read_delim(paste0("input/",dir,"/ancombc-formula/lfc_slice.csv"), delim=",")
qval <- read_delim(paste0("input/",dir,"/ancombc-formula/q_val_slice.csv"), delim=",")
data_aggregated <- merge(lfc, qval, by = "id", suffixes = c(".lfc", ".qval"))

# Step 4: Merge the taxonomic information with the LFC and q-value data
merged_data <- left_join(asv_taxonomy, data_aggregated, by = c("ASV_ID" = "id"))

# Step 5: Merge the tree's tip labels (ASVs) with the merged data
data_merged <- left_join(data.frame(id = tree$tip.label), merged_data, by = c("id" = "ASV_ID"))

# Step 6: Add a color column based on LFC (positive = red, negative = blue)
data_merged$lfc_color <- ifelse(data_merged$conditioncancer.lfc > 0, "red", "blue")

# Step 7: Add size based on q-value (larger for more significant)
data_merged$point_size <- ifelse(data_merged$conditioncancer.qval < 0.05, 0.5, 0.1)

# Step 8: Plot
p <- ggtree(tree, layout='circular', branch.length = "none")

p_data <- p$data %>%
  left_join(data_merged, by = c("label" = "id"))

p %+% p_data +
  geom_tiplab(size = 0.5) +
  geom_point(aes(color = lfc_color, size = point_size), show.legend = FALSE) +
  theme_tree()
