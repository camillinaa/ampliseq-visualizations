#!/usr/bin/env Rscript

# Phylogenetic trees
# Author: Camilla Callierotti - IU2 OMICS HT 
# Date: 03 March 2025

# ----- Load libraries -----
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(ggtree))



#----- Input files -----
parser <- OptionParser(usage = "%prog tax abundance tree meta contrast",
                       description = "-----------------------------------------------------------
  Description: This script performs diverity plots.
  Required input:
     1. Taxonomy file: ASV-level taxonomy file from DADA2 (dada2/ASV_tax.{database_name}.tsv). !!! CONSIDER USING /qiime2/abundance_tables/abs-abund-table-6.tsv
     2. Feature table: ASV-level feature abundance table from QIIME2 (qiime2/abundance_tables/feature-_table.tsv).
     2. Tree file: Newick tree file from QIIME2 (qiime2/phylogenetic_tree/tree.nwk).
     3. Metadata file: Metadata file with sample information.
     4. Contrast: is a string of format test:reference, with words separated by colon (:), where:
           - 'test' is the name of the factor level to be tested (e.g. treat). 
           - 'reference' is the name of the factor level to be used as reference (e.g. ctrl).
                        -----------------------------------------------------------")
arguments <- parse_args(parser, args <- commandArgs(trailingOnly = TRUE), positional_arguments = TRUE)
opt <- arguments$options

# Check if required arguments are provided
if (length(arguments$args) < 4) {
  stop("Error: taxonomy, tree, metadata and contrast files are required.")
}
# Use the arguments
dada2_asv_tax <- arguments$args[1]
abundace <- arguments$args[2]
tree_filename <- arguments$args[3]
meta_filename <- arguments$args[4]
contrast <- arguments$args[5]
print(paste("Contrast is:", contrast))

parts <- str_split(contrast, ":", simplify = TRUE)
formatted <- str_c(str_to_title(parts[1]), "vs", str_to_title(parts[2]), sep = " ")
formatted

#----- Prepare data -----
meta <- fread(meta)
feature <- fread(abundance)
asv_species <- fread(tax)

asv_species[asv_species == ""] <- NA
asv_species$taxonomy <- apply(asv_species[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")], 1, function(row) paste(na.omit(row), collapse = "; ")) # create taxonomy string
asv_species <- asv_species[asv_species$ASV_ID %in% tree$tip.label]
asv_species$Phylum[is.na(asv_species$Phylum)] <- "Unknown Phylum"

#----- Functions -----
# Update tip labels with taxonomy
tip_labels_to_taxonomy <- function(tree, asv_species) {
  tree$ASV <- tree$tip.label # Keep original ASVs
  tree$tip.label <- asv_species$taxonomy[match(tree$tip.label, asv_species$ASV_ID)] # Replace labels with taxonomy
  return(tree)
}
# Revert tip labels back to ASVs
tip_labels_to_asv <- function(tree) {
  if (!"ASV" %in% names(tree)) {
    stop("The tree does not contain ASV labels. Make sure to have run update_tip_labels first.")
  }
  tree$tip.label <- tree$ASV
  return(tree)
}

#----- Tree 1: Basic tree -----
# Change tip labels to ASV taxonomy
tree <- tip_labels_to_taxonomy(tree, asv_species)

# Colour tree by phylum
phyla_list <- list()
for (phylum in unique(asv_species$Phylum)) {
  taxonomies <- asv_species$taxonomy[asv_species$Phylum == phylum]
  phyla_list[[phylum]] <- taxonomies
}

# Create ggtree object
tree_coloured <- groupOTU(tree, phyla_list, 'Phylum')
t1 <- ggtree(tree_coloured, layout = 'daylight', branch.length = 'none', size = 1.5) +
  aes(color = Phylum) +
  theme(legend.position = "right", legend.text = element_text(size = 6))

# Save tree
ggsave("./nfdata/qiime2/tree_phyla.png", width = 100, height = 100, units = "cm", limitsize = FALSE)
ggsave("./nfdata/qiime2/tree_phyla.pdf", width = 100, height = 100, units = "cm", limitsize = FALSE)

#----- Tree 2: Circular tree with heatmap -----
# Create heatmap dataframe (from feature table using ASVs for join) !! Figure out how to use treatment variable from command line
feature_long <- melt(feature, id.vars = "#OTU ID", variable.name = "ID", value.name = "abundance") # reshape to long format
merged_feature <- merge(feature_long, meta[, .(ID, formatted[2])], by = "ID", all.x = TRUE) # merge with treatment
merged_feature_mean <- merged_feature[, .(mean_abund = mean(abundance)), by = .(`#OTU ID`, formatted[2])] # compute mean abundace for otu and treatment (not sample)
heatmap <- dcast(merged_feature_mean, `#OTU ID` ~ formatted[2], value.var = "mean_abund", fill = 0)
heatmap[`#OTU ID` %in% tree$ASV] # check that all OTUs in the heatmap are in the tree and vice versa
heatmap_matrix <- as.matrix(heatmap[, -1])
rownames(heatmap_matrix) <- heatmap$`#OTU ID`

# Change tip labels back to ASVs
tree <- tip_labels_to_asv(tree)

# Create ggtree object
circ <- ggtree(tree, layout = "circular", branch.length = 'none')
circ_heatmap <- gheatmap(circ, heatmap_matrix, offset=0, width=.2,
               colnames_angle=95, colnames_offset_y = .25) +
  scale_fill_gradient(low = "#C6D4F9", high = "#F1A7C2", name="Abundance")
circ_heatmap

# Save tree
ggsave("./nfdata/qiime2/tree_heatmap.png", plot = circ_heatmap, width = 100, height = 100, units = "cm", limitsize = FALSE)
ggsave("./nfdata/qiime2/tree_heatmap.pdf", plot = circ_heatmap, width = 100, height = 100, units = "cm", limitsize = FALSE)

##----- Tree 3: Tree with evolutionary distances -----
# Change tip labels to ASV taxonomy
tree <- tip_labels_to_taxonomy(tree, asv_species)

# Create ggtree object
horiz <- ggtree(tree) +
  theme_tree2() +
  geom_tiplab((label=tree$tip.label), size=3) +
  labs(caption="Evolutionary Distance")
horiz

# Save tree
ggsave("./nfdata/qiime2/tree_distances.png", plot = horiz, width = 100, height = 100, units = "cm", limitsize = FALSE)
ggsave("./nfdata/qiime2/tree_distances.pdf", plot = horiz, width = 100, height = 100, units = "cm", limitsize = FALSE)

# Create plotly object
#plotly::ggplotly(horiz, tooltip = c("label"))

# users can change branch length stored in tree object by using rescale_tree() function provided by the treeio package

#----- Tree 4: Tree with msa plot (TBD) -----
msaplot(p=ggtree(tree), fasta=paste0("data/",dir,"/",msa_fasta), window=c(150, 175))