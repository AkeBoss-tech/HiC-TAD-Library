#!/bin/bash
#SBATCH --job-name=hic_viz_sox11
#SBATCH --output=hic_viz_%j.out
#SBATCH --error=hic_viz_%j.err
#SBATCH --time=02:00:00            # You don't need 48 hours for this. 2 is plenty.
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4          # Cooler benefits from a few cores
#SBATCH --mem=64G                  # 64GB is safe for high-res matrices
#SBATCH --partition=main           # Use 'main' or 'mem' instead of 'gpu'
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ad2046@scarletmail.rutgers.edu

# ==========================================
# Amarel Cluster Hi-C Visualization Job
# ==========================================

echo "Job ID: $SLURM_JOB_ID"
echo "Working Directory: $(pwd)"

# --- 1. ENVIRONMENT SETUP ---
# You need to make sure you have 'cooler' and 'matplotlib' installed.
# If you haven't created a specific environment for this yet, 
# uncomment the lines below to create one on the fly (first time only).

source ~/.bashrc
# conda create -n hic_env python=3.9 -y
# conda activate hic_env
# pip install cooler matplotlib numpy

# Assuming you have an environment ready:
conda activate hic_env || echo "WARNING: Check your conda environment name"

# --- 2. RUN ANALYSIS ---
echo "Starting Visualization Script..."

# python visualize_tads.py
# Make sure the python script is in the same folder as this script
python visualize_tads.py

echo "Job Complete. Check the .png files in this directory."