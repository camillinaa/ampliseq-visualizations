#!/bin/bash
#SBATCH --job-name=tree-annotation   # Name of the job
#SBATCH --output=annotation-%j.out   # Standard output and error log
#SBATCH --error=annotation-%j.err    # Error log file
#SBATCH --time=01:00:00              # Maximum time for the job (HH:MM:SS)
#SBATCH --partition=cpuq             # Queue name (partition)
#SBATCH --ntasks=1                   # Number of tasks
#SBATCH --cpus-per-task=1            # Number of CPUs per task
#SBATCH --mem=4G                     # Memory allocation

# Activate the conda environment
source /facility/nfdata-omics/miniconda3/bin/activate /home/camilla.callierotti/r-ggtree

# Load R 
module load R

# Run the R script
Rscript /home/camilla.callierotti/microbiome/tree_annotation_script/annotation.R

# Deactivate conda environment
conda deactivate