"""
This script parses a KEGG (Kyoto Encyclopedia of Genes and Genomes) file to extract 
categories, modules, and entries, and saves the parsed data as a Pandas DataFrame.

It also parses Enzyme Commission (EC) numbers to extract enzyme classes and subclasses, 
and saves the parsed data as Pandas DataFrames.

Functions:
    parse_kegg_file(filepath):
        Reads a KEGG file and extracts hierarchical information about categories, 
        modules, and entries. Returns a Pandas DataFrame with columns:
        - "Entry": KEGG entry IDs (e.g., K numbers like K00001).
        - "Module": Module names associated with the entries.
        - "Category": High-level categories associated with the modules.

Usage:
    The script reads a KEGG file located at "../input/general/ko00001.keg", parses 
    its content, and saves the resulting DataFrame as a pickle file at 
    "../intermediate/kegg-categories.pkl".

    Additionally, it reads an enzyme classification file located at 
    "../input/general/enzclass.txt", parses its content to extract enzyme classes 
    and subclasses, and saves the resulting DataFrames as pickle files in 
    "../intermediate/ec-classes.pkl" and "../intermediate/ec-subclasses.pkl".
"""

import pandas as pd
import re


#--- Parse KEGG file for categories, modules, and entries ---


def parse_kegg_file(filepath):
    with open(filepath, "r") as file:
        kegg_data = file.readlines()
    
    entries = []
    modules = []
    categories = []
    
    module = None
    category = None
    
    for line in kegg_data:

        # Level A - Category (Starts with "A091")
        if line.startswith("A091"):
            category = line[1:].strip() # remove A 
            module = None # reset module when new category starts

        # Level B - Module
        elif line.startswith("B") and category:
            module = line.strip().split("  ", 1)[-1]

        # Level D - Entries
        elif line.startswith("D") and category and module:
            match = re.search(r"K\d{5,6}", line)
            if match:
                entries.append(match.group())
                modules.append(module)
                categories.append(category)
    
    return pd.DataFrame({"Entry": entries, "Module": modules, "Category": categories})

kegg_categories = parse_kegg_file("../input/general/ko00001.keg")
kegg_categories.to_pickle("../intermediate/kegg-categories.pkl")


#--- Parse Enzyme Commission Numbers (EC) ---


class_mapping = []
subclass_mapping = []

with open("input/general/enzclass.txt", "r") as f:
    for line in f:
        line = line.strip()

        # Match enzyme classes (e.g., "1. -.-  Oxidoreductases.")
        match_class = re.match(r"^(\d+)\.\s+-\.\s+-\.-\s+(.*)", line)

        # Match enzyme subclasses (e.g., "1. 1. -.-  Acting on the CH-OH group of donors.")
        match_subclass = re.match(r"^(\d+)\.\s*(\d+)\.\s+-\.-\s+(.*)", line)

        if match_class:
            class_mapping.append([match_class.group(1), match_class.group(2)])

        if match_subclass:
            subclass_code = f"{match_subclass.group(1)}.{match_subclass.group(2)}"
            subclass_mapping.append([subclass_code, match_subclass.group(3)])

class_df = pd.DataFrame(class_mapping, columns=["class", "class_name"])
subclass_df = pd.DataFrame(subclass_mapping, columns=["subclass", "subclass_name"])
class_df.to_pickle("intermediate/ec-classes.pkl")
subclass_df.to_pickle("intermediate/ec-subclasses.pkl")