
#
#   +----------------------------+
#   |  Load libraries & data     |
#   +----------------------------+
#   

message("Starting R script...")

# Import libraries 

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpubr)   # For statistical tests in ggplot
library(forcats)  # For ordering categories
library(ggpubr)   # For more customized layouts


# set vars 

dir = "pipeline_test"
ko_filename = "KO_pred_metagenome_unstrat_descrip.tsv"
ec_filename = "EC_pred_metagenome_unstrat_descrip.tsv"
metacyc_filename = "METACYC_path_abun_unstrat_descrip.tsv"
meta_filename = "Metadata.tsv"
condition_col = "treatment1"
results_path = "~/ampliseq_functional_analysis/plots/"

# Load data

meta <- read.delim(paste0("functional_data/",dir,"/",meta_filename), header = TRUE)
ko <- read.delim(paste0("functional_data/",dir,"/",ko_filename), header = TRUE)
ec <- read.delim(paste0("functional_data/",dir,"/",ec_filename), header = TRUE)
metacyc <- read.delim(paste0("functional_data/",dir,"/",metacyc_filename), row.names = 1, header = TRUE)

#
#   +----------------------------+
#   |  Prepare KO categories     |
#   +----------------------------+
#   

if (file.exists("ko_categories.rds")) {
  ko_categories <- readRDS("ko_categories.rds")
} else {
  kegg_data <- readLines(paste0("functional_data/",dir,"/ko00001.keg"))
  
  ko_terms <- c()
  categories <- c()
  
  for (i in 1:length(kegg_data)) {
    # Identify the category lines (e.g., "B  09101 Carbohydrate metabolism")
    if (grepl("^B", kegg_data[i])) {
      category <- str_trim(sub("^B\\s+", "", kegg_data[i]))
    }
    
    # Identify KO term lines (e.g., "D      K00844  HK; hexokinase [EC:2.7.1.1]")
    if (grepl("^D", kegg_data[i])) {
      ko_term <- str_trim(sub("^D\\s+", "", kegg_data[i]))
      ko_id <- str_extract(ko_term, "K\\d{5,6}")
      ko_terms <- c(ko_terms, ko_id)
      categories <- c(categories, category)
    }
  }
  
  ko_categories <- data.frame(KO = ko_terms, Category = categories, stringsAsFactors = FALSE)
  head(ko_categories)
  unique(ko_categories$Category)
  saveRDS(ko_categories, "ko_categories.rds")
}
  
#
#   +--------------------------+
#   |  Map Picrust results     |
#   +--------------------------+
#   

# Convert KEGG ontology abundances to long format and add metadata

ko_long <- ko %>%
  pivot_longer(cols = starts_with("sampleID"), 
               names_to = "Sample", 
               values_to = "Abundance") %>%
  left_join(ko_categories, by = c("function." = "KO"), 
            relationship = "many-to-many") %>%
  left_join(meta %>%
              select(ID, treatment1), by = c("Sample" = "ID"))

ko_long$treatment1 <- factor(ko_long$treatment1) # Convert to factor for statistical test in ggplot
ko_long$Group <- factor(ko_long$Group, levels = names(group_colors)) # Convert to factor for coloring in ggplot

# Group KO categories into groups

ko_long <- ko_long %>%
  mutate(
    Group = case_when(
      Category %in% c("09101 Carbohydrate metabolism", "09103 Lipid metabolism", "09105 Amino acid metabolism",
                      "09108 Metabolism of cofactors and vitamins", "09102 Energy metabolism",
                      "09109 Metabolism of terpenoids and polyketides", "09104 Nucleotide metabolism",
                      "09107 Glycan biosynthesis and metabolism", "09111 Xenobiotics biodegradation and metabolism",
                      "09110 Biosynthesis of other secondary metabolites") ~ "Metabolism",
      
      Category %in% c("09121 Transcription", "09122 Translation", "09124 Replication and repair", "09123 Folding, sorting and degradation") ~ "Genetic Information Processing",
      
      Category %in% c("09132 Signal transduction", "09141 Transport and catabolism", "09131 Membrane transport",
                      "09145 Cellular community - prokaryotes") ~ "Environmental Information Processing",
      
      Category %in% c("09142 Cell motility", "09143 Cell growth and death") ~ "Cellular Processes",
      
      Category %in% c("09151 Immune system", "09152 Endocrine system", "09154 Digestive system", "09156 Nervous system",
                      "09155 Excretory system", "09153 Circulatory system", "09158 Development and regeneration",
                      "09159 Environmental adaptation") ~ "Organismal Systems",
      
      Category %in% c("09161 Cancer: overview", "09162 Cancer: specific types", "09163 Immune disease",
                      "09164 Neurodegenerative disease", "09165 Substance dependence", "09166 Cardiovascular disease",
                      "09167 Endocrine and metabolic disease", "09171 Infectious disease: bacterial",
                      "09172 Infectious disease: viral", "09174 Infectious disease: parasitic",
                      "09175 Drug resistance: antimicrobial", "09176 Drug resistance: antineoplastic",
                      "09149 Aging") ~ "Human Diseases",
      
      Category %in% c("09181 Protein families: metabolism", "09182 Protein families: genetic information processing",
                      "09183 Protein families: signaling and cellular processes", "09191 Unclassified: metabolism",
                      "09192 Unclassified: genetic information processing", "09193 Unclassified: signaling and cellular processes",
                      "09194 Poorly characterized") ~ "Unclassified",
      
      TRUE ~ "Other"
    )
  )

# Define colors for each group

group_colors <- c(
  "Metabolism" = "#A6CEE3",
  "Genetic Information Processing" = "#1F78B4",
  "Environmental Information Processing" = "#B2DF8A",
  "Cellular Processes" = "#33A02C",
  "Organismal Systems" = "#FB9A99",
  "Human Diseases" = "#E31A1C",
  "Unclassified" = "#6A3D9A"
)


#
#   +------------------------------------------+
#   |    Plot Picrust KEGG Ontology results    |
#   +------------------------------------------+
#   

comparisons <- list(c("a", "b"))

# Boxplots of KEGG ontology abundance by treatment for each CATEGORY

ggplot(ko_long, aes(x = Abundance, y = treatment1, fill = treatment1)) +  
  geom_boxplot(outlier.size = 0.2) + 
  facet_grid(Group ~ Category, switch = "y") +  # Facet by Category, move labels to the left
  theme_minimal() +
  labs(title = "Relative Abundance by Treatment According to KEGG Ontology Category",
       x = "Abundance", 
       y = "Treatment",
       fill = "Treatment") +
  theme(
    axis.text.y = element_text(size = 8),
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    axis.title.y = element_blank(),   # Optional: Customize axis titles if needed
    strip.text.y.left = element_text(angle = 0),  # Rotate facet labels to be horizontal
    strip.placement = "outside"  # Move facet labels to the left
  ) +
  scale_y_discrete(labels = NULL)

library(grid)

grid.ls(grid.force())  

grid.gedit(gPath("GRID.stripGrob", "GRID.text"),  
           just = "left", x = unit(0, "npc"))


ggsave("ko_categories.png", width = 10, height = 8, dpi = 300)

# Boxplots of KEGG ontology abundance by treatment for each MACRO GROUP

ggplot(ko_long, aes(x = Abundance, y = treatment1, fill = treatment1)) +  
  geom_boxplot(outlier.size = 0.2) + 
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif") +  
  facet_grid(Group ~ ., switch = "y") +  # Facet by Group, move labels to the left
  theme_minimal() +
  labs(title = "Relative Abundance by Treatment According to KEGG Ontology Group",
       x = "Abundance", 
       y = "Treatment",
       fill = "Treatment") +
  theme(
    axis.text.y = element_text(size = 8),
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    axis.title.y = element_blank(),   # Optional: Customize axis titles if needed
    strip.text.y.left = element_text(angle = 0),  # Rotate facet labels to be horizontal
    strip.placement = "outside"  # Move facet labels to the left
  ) +
  scale_y_discrete(labels = NULL)

ggsave("ko_groups.png", width = 10, height = 8, dpi = 300)

# Functional analysis breakdown 


ggplot(ko_long, aes(x = Abundance, y = treatment1, fill = treatment1)) +  
  geom_boxplot(outlier.size = 0.2) + 
  facet_wrap(~ Group, strip.position = "left") +  # Facet by Group, move labels to the left
  theme_minimal() +
  labs(title = "Relative Abundance by Treatment According to KEGG Ontology Category",
       x = "Abundance", 
       y = "Treatment",
       fill = "Treatment") +
  theme(
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    axis.title.y = element_blank(),  # Remove y-axis title
    strip.text.y.left = element_text(size = 10, face = "bold", color = "white"),  # Set facet text color to white
    strip.background.x = element_rect(fill = group_colors[levels(factor(ko_long$Group))]),  # Set background color per group
    strip.placement = "outside"  # Move facet labels to the left
  ) +
  scale_y_discrete(labels = NULL)  # Remove the 'a' and 'b' labels on the y-axis

ggsave("functional_plot.png", p, width = 10, height = 12, dpi = 300)



##### trials


