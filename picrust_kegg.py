import os 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import re

# set vars 
cwd = os.getcwd()
dir = "pipeline-test"
kegg_filename = "KO_pred_metagenome_unstrat_descrip.tsv"
ec_filename = "EC_pred_metagenome_unstrat_descrip.tsv"
metacyc_filename = "METACYC_path_abun_unstrat_descrip.tsv"
meta_filename = "Metadata.tsv"
condition_col = "treatment1"

# Load Data
meta = pd.read_csv(f"{cwd}/input/{dir}/Metadata.tsv", sep="\t")
kegg = pd.read_csv(f"{cwd}/input/{dir}/picrust/KO_pred_metagenome_unstrat_descrip.tsv", sep="\t", header=0)

# Process KEGG ontology file to create kegg_categories
def parse_kegg_file(filepath):
    with open(filepath, "r") as file:
        kegg_data = file.readlines()
    
    kegg_terms = []
    categories = []
    category = None
    
    for line in kegg_data:
        if line.startswith("B"):
            category = line.strip().split("  ", 1)[-1]
        elif line.startswith("D") and category:
            match = re.search(r"K\d{5,6}", line)
            if match:
                kegg_terms.append(match.group())
                categories.append(category)
    
    return pd.DataFrame({"KO": kegg_terms, "Category": categories})

kegg_categories = parse_kegg_file("input/general/ko00001.keg")
kegg_categories.to_pickle("intermediate/kegg-categories.pkl")

# Convert kegg data to long format
value_columns = [col for col in kegg.columns if col not in ["function", "description"]]  # Get all columns that are not function or description
kegg_long = kegg.melt(id_vars=["function"], value_vars=value_columns, var_name="Sample", value_name="Abundance")  # Convert to long format
kegg_long = kegg_long.merge(kegg_categories, left_on="function", right_on="KO", how="left")  # Add KEGG category to dataframe
kegg_long = kegg_long.merge(meta[["ID", "treatment1"]], left_on="Sample", right_on="ID", how="left")  # Add treatment to dataframe
kegg_long["treatment1"] = kegg_long["treatment1"].astype("category")  # Convert treatment to categorical variable

# Assign groups based on categories
def assign_group(category):
    metabolism = ["09101 Carbohydrate metabolism", "09103 Lipid metabolism", "09105 Amino acid metabolism", 
                  "09108 Metabolism of cofactors and vitamins", "09102 Energy metabolism"]
    genetic_info = ["09121 Transcription", "09122 Translation", "09124 Replication and repair"]
    env_info = ["09132 Signal transduction", "09141 Transport and catabolism"]
    cell_processes = ["09142 Cell motility", "09143 Cell growth and death"]
    organismal = ["09151 Immune system", "09152 Endocrine system"]
    human_diseases = ["09161 Cancer: overview", "09162 Cancer: specific types"]
    unclassified = ["09181 Protein families: metabolism", "09182 Protein families: genetic information processing"]
    
    if category in metabolism:
        return "Metabolism"
    elif category in genetic_info:
        return "Genetic Information Processing"
    elif category in env_info:
        return "Environmental Information Processing"
    elif category in cell_processes:
        return "Cellular Processes"
    elif category in organismal:
        return "Organismal Systems"
    elif category in human_diseases:
        return "Human Diseases"
    elif category in unclassified:
        return "Unclassified"
    else:
        return "Other"

kegg_long["Group"] = kegg_long["Category"].apply(assign_group)
kegg_long["Category"] = kegg_long["Category"].fillna("").astype(str)

# Define colors for each group
group_colors = {
    "Metabolism": "#A6CEE3",
    "Genetic Information Processing": "#1F78B4",
    "Environmental Information Processing": "#B2DF8A",
    "Cellular Processes": "#33A02C",
    "Organismal Systems": "#FB9A99",
    "Human Diseases": "#E31A1C",
    "Unclassified": "#6A3D9A"
}

# Create figure
fig, ax = plt.subplots(figsize=(12, 8))

# Plot the boxplot (Horizontal boxplot)
sns.boxplot(x="Abundance", y="Category", hue="treatment1", data=kegg_long, order=kegg_long["Category"].unique(), ax=ax, showfliers=False, palette=["#E41A1C", "#377EB8"])

# Rotate x labels (No need to rotate as it's horizontal now)
plt.xticks(rotation=45, ha='right')


# Dynamically generate group bounds based on the categories in kegg_long
group_bounds = []
group_labels = []

# Sort categories based on their group
for group, group_data in kegg_long.groupby("Group"):
    group_labels.append(group)
    
    # Get sorted list of categories for the current group
    sorted_categories = sorted(group_data["Category"].unique())
    start_idx = len(group_bounds)  # Start index for the current group
    end_idx = start_idx + len(sorted_categories)  # End index for the current group
    
    # Add the bounds for the group
    group_bounds.append((start_idx, end_idx))

# Add background colors for each group (shade horizontal regions along y-axis)
for i, (group, (y_start, y_end), color) in enumerate(zip(group_labels, group_bounds, group_colors.values())):
    ax.axhspan(y_start - 0.5, y_end - 0.5, color=color, alpha=0.3, zorder=0)  # Background color
    ax.text(max(kegg_long['Abundance']) + 0.1, (y_start + y_end) / 2 - 0.5, group, 
            ha='left', va='center', fontsize=12, fontweight='bold')

# Title and labels
plt.title('KEGG Ontology Relative Abundances by Treatment', fontsize=14)
plt.xlabel('Pathway Abundance Value')
plt.ylabel('KEGG Category')

# Adjust layout
plt.tight_layout()

# Save the plot
plt.show()
#plt.savefig("results/picrust/kegg.png", dpi=300)

# Close to free up memory
plt.close()
# %%
