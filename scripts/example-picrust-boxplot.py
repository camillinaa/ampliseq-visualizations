import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

# Example data structure
data = pd.DataFrame({
    'category': np.random.choice([
        'D-Glutamine and D-glutamate metabolism', 'D-Alanine metabolism', 
        'Pentose phosphate pathway', 'Fructose and mannose metabolism',
        'Carbon fixation in photosynthetic organisms', 'Sulfur metabolism'
    ], 100),
    'abundance': np.random.rand(100) * 2,
    'treatment': np.random.choice(['Control', 'GMCA'], 100),
    'group': np.random.choice([
        'Amino Acid Metabolism', 'Carbohydrate and Lipid Metabolism', 'Energy Metabolism'
    ], 100)
})

# Define order of categories
category_order = [
    'D-Glutamine and D-glutamate metabolism', 'D-Alanine metabolism',
    'Pentose phosphate pathway', 'Fructose and mannose metabolism',
    'Carbon fixation in photosynthetic organisms', 'Sulfur metabolism'
]

# Define groups and their corresponding background colors
group_labels = ['Amino Acid Metabolism', 'Carbohydrate and Lipid Metabolism', 'Energy Metabolism']
group_bounds = [(0, 2), (2, 4), (4, 6)]  # Define start and end x-axis positions
group_colors = ['#A6CEE3', '#FDBF6F', '#B2DF8A']  # Soft pastel colors

# Create figure
fig, ax = plt.subplots(figsize=(10, 6))

# Plot the boxplot
sns.boxplot(x='category', y='abundance', hue='treatment', data=data,
            order=category_order, ax=ax, showfliers=False)

# Rotate x labels
plt.xticks(rotation=45, ha='right')

# Add background colors for groups
for i, (group, (x_start, x_end), color) in enumerate(zip(group_labels, group_bounds, group_colors)):
    ax.axvspan(x_start - 0.5, x_end - 0.5, color=color, alpha=0.3, zorder=0)  # Background color
    ax.text((x_start + x_end) / 2 - 0.5, max(data['abundance']) + 0.1, group, 
            ha='center', va='bottom', fontsize=12, fontweight='bold')

# Title and labels
plt.title('Pathway Abundance for Different Treatments', fontsize=14)
plt.xlabel('KEGG Category')
plt.ylabel('Pathway Abundance Value')

# Adjust layout
plt.tight_layout()

# Save the plot
plt.savefig('example_picrust_boxplot.png', dpi=300)

# Close to free up memory
plt.close()
