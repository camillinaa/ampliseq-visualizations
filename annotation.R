# 
#        _           _                             _   _        _                 
#       | |         | |                           | | (_)      | |                
#  _ __ | |__  _   _| | ___   __ _  ___ _ __   ___| |_ _  ___  | |_ _ __ ___  ___ 
# | '_ \| '_ \| | | | |/ _ \ / _` |/ _ \ '_ \ / _ \ __| |/ __| | __| '__/ _ \/ _ \
# | |_) | | | | |_| | | (_) | (_| |  __/ | | |  __/ |_| | (__  | |_| | |  __/  __/
# | .__/|_| |_|\__, |_|\___/ \__, |\___|_| |_|\___|\__|_|\___|  \__|_|  \___|\___|
# | |           __/ |         __/ |                                               
# |_|          |___/         |___/
# 
#   
# nf-core/ampliseq phylogenetic tree generation
# february 2025
# sources:
# https://yulab-smu.top/treedata-book/chapter2.html
# https://yulab-smu.top/treedata-book/related-tools.html#plotly
# 

#
#   +----------------------------+
#   |  Load libraries & data     |
#   +----------------------------+
#   

# Import libraries 

library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(ggtree)

# Set vars

dir = "pipeline_test"
condition_col = "treatment1"

# Load Newick tree

tree <- read.tree(paste0("data/",dir,"/tree.nwk"))

# Load taxonomy_qza

taxonomy_qza <- qiime2R::read_qza(paste0("data/",dir,"/taxonomy.qza"))$data
max_levels <- max(stringr::str_count(taxonomy_qza$Taxon, ";")) + 1
taxonomy <- taxonomy_qza %>%
  rename(ASV = Feature.ID, Taxonomy = Taxon) %>%
  select(ASV, Taxonomy) %>%
  tidyr::separate(Taxonomy, into = paste0("Rank", 1:max_levels), sep = ";", fill = "right")

# Load and prepare metadata

meta <- fread(paste0("data/",dir,"/Metadata.tsv"))
feature <- fread(paste0("data/",dir,"/feature-table.tsv"))
asv_species <- fread(paste0("data/",dir,"/ASV_tax.gtdb_R07-RS207.tsv"))

asv_species[asv_species == ""] <- NA
asv_species$taxonomy <- apply(asv_species[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")], 1, function(row) paste(na.omit(row), collapse = "; ")) # create taxonomy string
asv_species <- asv_species[asv_species$ASV_ID %in% tree$tip.label]
asv_species$Phylum[is.na(asv_species$Phylum)] <- "Unknown Phylum"

#
#   +----------------------------+
#   |  TREE 1 - Basic tree       |
#   +----------------------------+
#   

# Change tip labels to ASV taxonomy

tree$ASV <- tree$tip.label # keep ASVs in different object
tree$tip.label <- asv_species$taxonomy[match(tree$tip.label, asv_species$ASV_ID)] # change tip label to taxonomy

# Colour tree by phylum

phyla_list <- list()
for (phylum in unique(asv_species$Phylum)) {
  taxonomies <- asv_species$taxonomy[asv_species$Phylum == phylum]
  phyla_list[[phylum]] <- taxonomies
}

# Create ggtree object

t1 <- ggtree(tree, layout = 'daylight', branch.length = 'none', size=1.5) + geom_tiplab(size = 6)
t1_coloured <- groupOTU(t1, phyla_list, 'Phylum') +
  aes(color = Phylum) +
  theme(legend.position = "right", legend.text = element_text(size = 6))
t1_coloured

# Save tree

ggsave("trees/tree_phyla.pdf", width = 100, height = 100, units = "cm", limitsize = FALSE)

#
#   +--------------------------------------------+
#   |    TREE 2 - Circular tree with heatmap     |
#   +--------------------------------------------+
#   

# Create heatmap dataframe (from feature table using ASVs for join)

feature_long <- melt(feature, id.vars = "#OTU ID", variable.name = "ID", value.name = "abundance") # reshape to long format
merged_feature <- merge(feature_long, meta[, .(ID, treatment1)], by = "ID", all.x = TRUE) # merge with treatment
heatmap <- merged_feature[, .(mean_abund = mean(abundance, na.rm = TRUE)), by = c("#OTU ID", "treatment1")] # compute mean for otu/treatment pair
heatmap <- dcast(heatmap, `#OTU ID` ~ treatment1, value.var = "mean_abund", fill = 0) # reshape back to wide format

heatmap[`#OTU ID` %in% tree$ASV] # check that all OTUs in the heatmap are in the tree and vice versa

heatmap[, 2:3] <- lapply(heatmap[, 2:3], as.numeric)
heatmap_matrix <- as.matrix(heatmap[, -1])
rownames(heatmap_matrix) <- heatmap$`#OTU ID`

# Change tip labels back to ASVs

tree$tip.label <- tree$ASV

# Create ggtree object

circ <- ggtree(tree, layout = "circular", branch.length = 'none')
circ_heatmap <- gheatmap(circ, heatmap_matrix, offset=0, width=.2,
               colnames_angle=95, colnames_offset_y = .25) +
  scale_fill_gradient(low = "#C6D4F9", high = "#F1A7C2", name="Abundance")
circ_heatmap

# Save tree

ggsave("trees/tree_heatmap.pdf", width = 100, height = 100, units = "cm", limitsize = FALSE)

#
#   +--------------------------------------------+
#   |  TREE 3 - Tree with evolutionary distances |
#   +--------------------------------------------+
#   

horiz <- ggtree(tree) + 
  theme_tree2() + 
  geom_tiplab(size=3) + 
  labs(caption="Evolutionary Distance")
  
# users can change branch length stored in tree object by using rescale_tree() function provided by the treeio package

