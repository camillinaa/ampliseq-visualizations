# Import libraries 

library(dplyr)
library(tidyr)
library(ggtree)
library(ggnewscale) # for multiple heatmaps

# Set vars

dir = "pipeline_test"
condition_col = "treatment1"

# Load Newick tree and taxonomy

tree <- read.tree(paste0("data/",dir,"/tree.nwk"))
taxonomy_qza <- qiime2R::read_qza(paste0("data/",dir,"/taxonomy.qza"))$data
max_levels <- max(stringr::str_count(taxonomy_qza$Taxon, ";")) + 1
taxonomy <- taxonomy_qza %>%
  rename(ASV = Feature.ID, Taxonomy = Taxon) %>%
  select(ASV, Taxonomy) %>%
  tidyr::separate(Taxonomy, into = paste0("Rank", 1:max_levels), sep = ";", fill = "right")

# Load metadata

metadata <- read.csv(paste0("data/",dir,"/level-9.csv"))

# Insert taxonomy in tree tip labels and clean

tree$ASV <- tree$tip.label
tree$tip.label <- taxonomy_qza$Taxon[match(tree$tip.label, taxonomy_qza$Feature.ID)]
tree$tip.label <- lapply(tree$tip.label, sub, pattern = ";;.*", replacement = "")



### TREE 1: basic phylum-coloured tree



# Colour tree by phylum

asv_phylum <- setNames(taxonomy$Rank2, taxonomy$ASV)
asv_phylum_tree <- asv_phylum[tree$ASV]
asv_phylum_tree[asv_phylum_tree == ""] <- "Unknown Phyla"
grp <- lapply(split(tree$tip.label, asv_phylum_tree), unlist, use.names = FALSE)

# Create ggtree object

daylight_tree <- ggtree(tree, layout = 'daylight', branch.length = 'none', size=1.5) + geom_tiplab(size = 6)
daylight_tree_coloured <- groupOTU(daylight_tree, grp, 'Phylum') + aes(color = Phylum) +
  theme(legend.position = "right", legend.text=element_text(size=6))
daylight_tree_coloured

# Save tree

ggsave("trees/circular_tree.pdf", width = 100, height = 100, units = "cm", limitsize = FALSE)



### TREE 2: Circular tree with heatmap of condition variable (no colored branches)


# Create heatmap dataframe

df <- meta
df_summary <- df %>%
  group_by(treatment1) %>%
  summarise(across(starts_with("Bacteria"), mean, na.rm = TRUE))
df <- t(df_summary)
colnames(df) <- df[1, ]
df <- df[-1, ]

rownames(df) <- gsub("\\.", ";", rownames(df))
rownames(df) <- sub(";;.*", "", rownames(df))



# Create ggtree object

circ <- ggtree(tree, layout = "circular")

circ_heatmap <- gheatmap(circ, df, offset=.8, width=.2,
               colnames_angle=95, colnames_offset_y = .25) +
  scale_fill_viridis_d(option="D", name="Condition")


### TREE 3: Horizontal tree with annotations and evolutionary distances

horiz <- ggtree(tree) + 
  geom_tiplab(size=3) + 
  theme_tree2() #+ 
  labs(caption="Evolutionary Distance")
  
# users can change branch length stored in tree object by using rescale_tree() function provided by the treeio package

