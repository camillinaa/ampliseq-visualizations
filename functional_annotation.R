
#
#   +----------------------------+
#   |  Load libraries & data     |
#   +----------------------------+
#   

message("Starting R script...")

# Import libraries 

library(dplyr)
library(stringr)
library(ggplot2)

# set vaggplot2# set vars (ampliseq test profile)

dir = "pipeline_test"
ko_filename = "KO_pred_metagenome_unstrat_descrip.tsv"
ec_filename = "EC_pred_metagenome_unstrat_descrip.tsv"
metacyc_filename = "METACYC_path_abun_unstrat_descrip.tsv"
meta_filename = "Metadata.tsv"
condition_col = "treatment1"
results_path = "~/ampliseq_functional_analysis/plots/"

# Load data

meta <- read.delim(paste0("data/",dir,"/",meta_filename), header = TRUE)
ko <- read.delim(paste0("data/",dir,"/",ko_filename), header = TRUE)
ec <- read.delim(paste0("data/",dir,"/",ec_filename), header = TRUE)
metacyc <- read.delim(paste0("data/",dir,"/",metacyc_filename), row.names = 1, header = TRUE)

#
#   +----------------------------+
#   |  Prepare KO categories     |
#   +----------------------------+
#   

kegg_data <- readLines(paste0("data/",dir,"/ko00001.keg"))

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

ko_data <- data.frame(KO = ko_terms, Category = categories, stringsAsFactors = FALSE)

head(ko_data)
unique(ko_data$Category)

#
#   +--------------------------+
#   |  Map Picrust results     |
#   +--------------------------+
#   

# Convert ko data to long format

ko_long <- ko %>%
  pivot_longer(cols = starts_with("sampleID"), 
               names_to = "Sample", 
               values_to = "Abundance") %>%
  left_join(ko_data, by = c("function." = "KO"), 
            relationship = "many-to-many") %>%
  left_join(meta %>%
              select(ID, treatment1), by = c("Sample" = "ID"))

#
#   +----------------------------+
#   |    Plot Picrust results    |
#   +----------------------------+
#   

# Boxplots of KO abundance by treatment for each KO category

ggplot(ko_long, aes(x = treatment1, y = Abundance, fill = treatment1)) +
  geom_boxplot() +
  facet_wrap(~ Category, scales = "free_y") +  # Use free_y to adjust y-axis scale per category
  theme_minimal() +
  labs(title = "KO Abundance by Treatment and Category",
       x = "Treatment Group",
       y = "Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

