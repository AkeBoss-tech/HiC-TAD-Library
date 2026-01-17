import cooler
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import os

# got data from 
# https://data.4dnucleome.org/experiment-set-replicates/4DNES14CNC1I/#details
# https://data.4dnucleome.org/files-processed/4DNFINNZDDXV/

# --- MAC CONFIGURATION ---
# Update this to point to where you downloaded the file
FILE_PATH = "/Users/akashdubey/Documents/CodingProjects/HiC-TAD-Library/data/raw/mouse_microc.mcool" 

# Micro-C has very high resolution. 
# 5000 (5kb) is great for seeing TADs. 10000 (10kb) is safer for memory.
RESOLUTION = 5000  

# Coordinates for Mouse (mm10 assembly)
# We add a buffer to see the boundaries around the genes
regions = {
    "Sox11_Chr12": "chr12:26,000,000-28,000,000", 
    "Mir9-2_Chr13": "chr13:83,500,000-84,500,000"
}

def plot_region(file_path, region_name, coordinates):
    print(f"Processing {region_name}...")
    
    # 1. Access the specific resolution inside the mcool file
    uri = f"{file_path}::resolutions/{RESOLUTION}"
    
    try:
        clr = cooler.Cooler(uri)
    except Exception as e:
        print(f"Error: Could not load file. Check the path! {e}")
        return

    # 2. Fetch the matrix (balanced/normalized)
    # The 'fetch' command only loads this tiny slice into RAM
    matrix = clr.matrix(balance=True).fetch(coordinates)

    # 3. Visualize
    plt.figure(figsize=(8, 8))
    
    # Use a log scale because the diagonal is heavily dominant
    im = plt.imshow(
        matrix, 
        cmap='RdYlBu_r', 
        interpolation='none',
        norm=colors.LogNorm(vmin=0.001, vmax=0.05) 
    )
    
    plt.title(f"{region_name}\n{coordinates}")
    plt.colorbar(im, label="Contact Frequency")
    
    # 4. Save
    output_filename = f"{region_name}_heatmap.png"
    plt.savefig(output_filename, dpi=150)
    print(f"Saved to {output_filename}")
    plt.close()

if __name__ == "__main__":
    if not os.path.exists(FILE_PATH):
        print(f"STOP: I can't find the file at {FILE_PATH}")
    else:
        for name, coords in regions.items():
            plot_region(FILE_PATH, name, coords)