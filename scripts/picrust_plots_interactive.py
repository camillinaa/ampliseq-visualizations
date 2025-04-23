"""
08_picrust_plots.py
===================
This script performs functional analysis and visualization of PiCRUST2 output data, 
including KEGG Orthology and Enzyme Commission analyses. It processes input files, 
generates long-format dataframes, and creates boxplots to visualize relative abundances
by treatment or condition.

Functions:
-----------
- parse_args(): Parses command-line arguments for input file paths and metadata details.
- main(): Main function to execute the script logic.
- parse_kegg_file(filepath): Processes KEGG orthology file to extract categories, modules, and entries.

Usage:
------
Run the script from the command line with the required arguments:
    python 08_picrust_plots.py --kegg <kegg_file> --ec <ec_file> --meta <metadata_file> --condition_col <condition_column>

Arguments:
----------
- --kegg: KEGG functional orthologs file (e.g., "KO_pred_metagenome_unstrat_descrip.tsv").
- --ec: Enzyme Commission file (e.g., "EC_pred_metagenome_unstrat_descrip.tsv").
- --meta: Metadata file containing sample information (e.g., "Metadata.tsv").
- --condition_col: Column name in the metadata file representing the condition or treatment (e.g., "treatment1").

Dependencies:
-------------
- pandas
- numpy
- matplotlib
- seaborn
- argparse
- re
"""

import argparse
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from io import StringIO


def parse_args():
    parser = argparse.ArgumentParser(description="Run PiCRUST2 Functional Analysis")
    parser.add_argument("--dir", required=True, help="Directory containing input files")
    parser.add_argument("--kegg", required=True, help="KEGG functional orthologs file")
    parser.add_argument("--ec", required=True, help="Enzyme Commission file")
    # parser.add_argument("--kegg_categories_file", required=True, help="Path to KEGG categories file") 
    parser.add_argument("--meta", required=True, help="Metadata file")
    parser.add_argument("--condition_col", required=True, help="Column name for condition/treatment in metadata")
    return parser.parse_args()


def plot_boxplot(data, x_col, y_col, hue_col, group_col, group_bounds, group_colors, group_labels, module_order, title, xlabel, ylabel, output_file):
    """
    Generalized function to create a boxplot with background shading for groups.

    Parameters:
    - data: DataFrame containing the data to plot.
    - x_col: Column name for x-axis values.
    - y_col: Column name for y-axis values.
    - hue_col: Column name for hue (e.g., treatment/condition).
    - group_col: Column name for grouping (e.g., categories or classes).
    - group_bounds: List of tuples indicating start and end bounds for each group.
    - group_colors: List of colors for each group.
    - title: Plot title.
    - xlabel: Label for x-axis.
    - ylabel: Label for y-axis.
    - output_file: File path to save the plot.
    """
    
    fig, ax = plt.subplots(figsize=(14, 16))

    # Remove whitespace in plot
    category_positions = range(len(data[y_col].unique()))
    ax.set_ylim(min(category_positions), max(category_positions))

    # Seaborn boxplot
    sns.boxplot(
        data=data,
        x=x_col,
        y=y_col,
        order=module_order,  
        hue=hue_col,
        flierprops=dict(marker='o', markerfacecolor='black', markersize=0.75),
        ax=ax
    )

    # Add background colors for each group
    for i, (group, (y_start, y_end), color) in enumerate(zip(data[group_col].unique(), group_bounds, group_colors)):
        print(f"Group: {group}, Bounds: ({y_start}, {y_end}), Color: {color}")
        ax.axhspan(y_start - 0.5, y_end - 0.5, color=color, alpha=0.3, zorder=0)  # Background color
        ax.text(ax.get_xlim()[1], (y_start + y_end) / 2 - 0.5, group,
                ha='left', va='center', fontsize=12, fontweight='bold')

    # Title and labels
    plt.title(title, fontsize=16, fontweight='bold')
    plt.xlabel(xlabel, fontweight='bold')
    plt.ylabel(ylabel, fontweight='bold')
    plt.legend(title=args.condition_col.capitalize(), loc='lower right')

    # Save plot
    plt.savefig(output_file)
    plt.show()


def process_kegg_data(kegg_file, meta_file, kegg_categories_pkl, condition_col):

    """
    Process KEGG data for plotting.

    Parameters:
    - kegg_file: Path to KEGG data file.
    - meta_file: Path to metadata file.
    - kegg_categories_pkl: Path to KEGG categories pickle file parsed by ontology_parser.py.
    - condition_col: Column name for condition/treatment in metadata.

    Returns:
    - kegg_long (pd.DataFrame): Processed DataFrame in long format, ready for plotting.
    - group_bounds (list of tuples): List of tuples indicating the start and end indices for each group in the y-axis.
    - group_colors (list of str): List of colors corresponding to each group for plotting.
    - group_labels (list of str): List of unique group labels (categories) for the y-axis.
    - module_order (list of str): Ordered list of KEGG modules for plotting.
    """
    
    meta = pd.read_csv(meta_file, sep="\t")
    kegg = pd.read_csv(kegg_file, sep="\t", header=0)
    kegg_categories = pd.read_pickle(kegg_categories_pkl)

    # Convert KEGG data to long format
    value_columns = [col for col in kegg.columns if col not in ["function", "description"]] # extract columns with abundance values
    kegg_long = kegg.melt(id_vars=["function"], value_vars=value_columns, var_name="Sample", value_name="Abundance") # melt to long format to get samples on a column
    kegg_long = kegg_long.merge(kegg_categories, left_on="function", right_on="Entry", how="left") # add description
    kegg_long = kegg_long.merge(meta[["ID", condition_col]], left_on="Sample", right_on="ID", how="left").drop(columns=["ID"]) # add metadata for sample
    kegg_long[condition_col] = kegg_long[condition_col].astype("category")
    kegg_long = kegg_long.dropna() # Remove NA values (to be consistent with group bounds and colors)

    # Compute bounds and colors
    df = kegg_long[["Category", "Module"]].drop_duplicates()
    df["Category"] = df["Category"].astype("object").fillna("Unassigned")
    group_labels = df.loc[df["Category"] != "Unassigned", "Category"].unique().tolist() # don't include unassigned hereon after
    group_counts = df["Category"].value_counts(sort=False).loc[group_labels].tolist()
    cumulative_counts = [0] + list(np.cumsum(group_counts).astype(int))
    group_bounds = [(cumulative_counts[i], cumulative_counts[i + 1]) for i in range(len(group_counts))]
    group_colors = ['#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99', '#E31A1C', '#FDBF6F', '#FF7F00'][:len(group_labels)]

    # Create y-axis order: KEGG modules grouped by Category
    kegg_sorted = kegg_long[kegg_long["Category"] != "Unassigned"].copy()
    kegg_sorted["Category"] = pd.Categorical(kegg_sorted["Category"], categories=group_labels, ordered=True)
    kegg_sorted = kegg_sorted.sort_values(["Category", "Module"])
    module_order = kegg_sorted["Module"].unique().tolist()

    # Debugging: Print intermediate variables
    print("DataFrame (df):")
    print(df) # CORRECT

    print("\nGroup Labels:")
    print(group_labels)

    print("\nGroup Counts:")
    print(group_counts)

    print("\nGroup Bounds:")
    print(group_bounds)

    return kegg_long, group_bounds, group_colors, group_labels, module_order


def process_ec_data(ec_file, meta_file, class_pkl, subclass_pkl, condition_col):
    """
    Process EC data for plotting.

    Parameters:
    - ec_file: Path to EC data file.
    - meta_file: Path to metadata file.
    - class_pkl: Path to enzyme class mapping pickle file parsed by ontology_parser.py.
    - subclass_pkl: Path to enzyme subclass mapping pickle file parsed by ontology_parser.py.
    - condition_col: Column name for condition/treatment in metadata.

    Returns:
    - Processed DataFrame ready for plotting.
    - Group bounds and colors for plotting.
    - Group labels (enzyme class names).
    - Module order (subclass_name list) for y-axis sorting.
    """
    
    meta = pd.read_csv(meta_file, sep="\t")
    ec = pd.read_csv(ec_file, sep="\t", header=0)
    class_df = pd.read_pickle(class_pkl)
    subclass_df = pd.read_pickle(subclass_pkl)

    # Add class and subclass to ec
    ec["class"] = ec["function"].str.extract(r"EC:(\d)\.\d")
    ec["subclass"] = ec["function"].str.extract(r"EC:(\d\.\d)\..")
    ec = ec.merge(class_df, how="left", on="class")
    ec = ec.merge(subclass_df, how="left", on="subclass")

    # Summarize by subclass
    numeric_cols = ec.columns[2:-4] # don't include the first two because they are function and description, don't include the last four because they are class and subclass numbers and names
    ec_summarized = ec.groupby(["class", "class_name", "subclass", "subclass_name"])[numeric_cols].sum().reset_index() # keep only up to subclass level, collapsing function description

    # Melt the DataFrame to long format
    ec_long = ec_summarized.melt(
        id_vars=["class", "class_name", "subclass", "subclass_name"],
        var_name="sample",
        value_name="abundance"
    )

    # Merge metadata
    ec_long = ec_long.merge(meta[["ID", condition_col]], how="left", left_on="sample", right_on="ID").drop("ID", axis=1)
    
    # Clean dots in class and subclass names
    for colname in ["class_name", "subclass_name"]:
        ec_long[colname] = ec_long[colname].str.replace(".", "")

    # Compute group bounds and colors
    unique_groups = ec_long[["class_name", "subclass_name"]].drop_duplicates().reset_index()

    # Ensure the 'class_name' and 'subclass_name' are treated as categorical variables and ordered correctly.
    ec_sorted = ec_long.copy()

    # Step 1: Assign 'class_name' and 'subclass_name' as categorical types to preserve the order
    ec_sorted['class_name'] = pd.Categorical(ec_sorted['class_name'], categories=ec_sorted['class_name'].unique(), ordered=True)
    ec_sorted['subclass_name'] = pd.Categorical(ec_sorted['subclass_name'], categories=ec_sorted['subclass_name'].unique(), ordered=True)

    # Step 2: Sort values by 'class_name' first, then by 'subclass_name' to preserve hierarchical structure
    ec_sorted = ec_sorted.sort_values(by=["class_name", "subclass_name"], ascending=[True, True])

    # Step 3: Correctly assign 'module_order' based on the sorted order of 'subclass_name'
    module_order = ec_sorted['subclass_name'].unique().tolist()

    print("\nY-axis (subclass) module_order:")
    print(module_order)

    # Step 4: Create a mapping from subclass to class (which is already correct from your 'subclass_to_class' dataframe)
    subclass_to_class = ec_sorted.drop_duplicates(subset=["subclass_name"])[["subclass_name", "class_name"]]

    print("\nSubclass to Class mapping:")
    print(subclass_to_class)

    # Step 5: Explicitly order the 'class_name' for consistent grouping in the plot
    unique_groups = ec_sorted[["class_name", "subclass_name"]].drop_duplicates().reset_index()

    # Make sure the 'class_name' is ordered based on the EC hierarchy
    module_order_classes = unique_groups["class_name"].unique().tolist()

    # Step 6: Set the 'class_name' as categorical with the explicit order from 'module_order_classes'
    unique_groups['class_name'] = pd.Categorical(unique_groups['class_name'], categories=module_order_classes, ordered=True)

    # Step 7: Now group by 'class_name' and compute the bounds for each group
    group_bounds = []
    group_colors = []
    group_labels = []

    color_map = {  
        "Oxidoreductases": "#A6CEE3",
        "Transferases": "#1F78B4",
        "Hydrolases": "#B2DF8A",
        "Lyases": "#33A02C",
        "Isomerases": "#FB9A99",
        "Ligases": "#E31A1C",
        "Translocases": "#FDBF6F"
    }

    # Compute the bounds and colors for each group
    for class_name, group_df in unique_groups.groupby("class_name"):
        start_idx = group_df["index"].min()  # Find the starting index of the group
        end_idx = group_df["index"].max() + 1  # Find the ending index, adjust for inclusive range
        group_bounds.append((start_idx, end_idx))
        group_colors.append(color_map[class_name])  # Use a predefined color map for the groups
        group_labels.append(class_name)

    # Step 8: Output the final bounds and labels to verify the ordering
    for label, bounds in zip(group_labels, group_bounds):
        print(f"Group: {label}, Bounds: {bounds}")

    # Step 9: Return the necessary values for further analysis or plotting
    print(group_bounds, group_labels)
    return ec_long, group_bounds, group_colors, group_labels, module_order



def main():
    args = parse_args()

    # Set file names from command-line arguments
    dir = args.dir
    kegg_filename = args.kegg
    ec_filename = args.ec
    meta_filename = args.meta
    kegg_categories_filename = args.kegg_categories_file
    condition_col = args.condition_col


if __name__ == "__main__":
    args = parse_args()

    # Process KEGG data
    kegg_long, kegg_group_bounds, kegg_group_colors, kegg_group_labels, kegg_module_order= process_kegg_data(
        kegg_file=args.kegg,
        meta_file=args.meta,
        kegg_categories_pkl="intermediate/kegg-categories.pkl",
        condition_col=args.condition_col
    )

    # Plot KEGG data
    plot_boxplot(
        data=kegg_long,
        x_col="Abundance",
        y_col="Module",
        hue_col=args.condition_col,
        group_col="Category",
        group_bounds=kegg_group_bounds,
        group_colors=kegg_group_colors,
        group_labels=kegg_group_labels,
        module_order=kegg_module_order,
        title=f"KEGG Orthology Relative Abundances by {args.condition_col.capitalize()}",
        xlabel="Pathway Abundance Value",
        ylabel="KEGG Orthology Category",
        output_file="results/picrust/picrust_kegg.pdf"
    )

    # Process EC data
    ec_long, ec_group_bounds, ec_group_colors, ec_group_labels, ec_module_order = process_ec_data(
        ec_file=args.ec,
        meta_file=args.meta,
        class_pkl="intermediate/ec-classes.pkl",
        subclass_pkl="intermediate/ec-subclasses.pkl",
        condition_col=args.condition_col
    )

    # Plot EC data
    plot_boxplot(
        data=ec_long,
        x_col="abundance",
        y_col="subclass_name",
        hue_col=args.condition_col,
        group_col="class_name",
        group_bounds=ec_group_bounds,
        group_colors=ec_group_colors,
        group_labels=ec_group_labels,
        module_order=ec_module_order,
        title=f"Enzyme Commission Relative Abundances by {args.condition_col.capitalize()}",
        xlabel="Abundance Value",
        ylabel="Enzyme Commission Subclass",
        output_file="results/picrust/picrust_ec.pdf"
    )
