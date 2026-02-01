import cooler
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import os
import bioframe
import cooltools
from matplotlib.gridspec import GridSpec

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
    "Mir9-2_Chr13": "chr13:83,500,000-84,500,000",
    "Compartments_Chr2": "chr2:0-50,000,000"
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
    output_filename = f"media/{region_name}_heatmap.png"
    plt.savefig(output_filename, dpi=150, bbox_inches='tight')
    print(f"Saved to {output_filename}")
    plt.close()

def plot_triangular_region(file_path, region_name, coordinates):
    print(f"Processing triangular {region_name}...")
    
    uri = f"{file_path}::resolutions/{RESOLUTION}"
    try:
        clr = cooler.Cooler(uri)
    except Exception as e:
        print(f"Error: Could not load file. {e}")
        return

    # 1. Fetch matrix
    matrix = clr.matrix(balance=True).fetch(coordinates)
    n = matrix.shape[0]

    # 2. Create coordinates for pcolormesh
    # Standard Hi-C triangle transformation:
    # x = (i + j) / 2
    # y = (j - i) / 2
    # We use n+1 to define the corners of the bins
    x = np.arange(n + 1)
    y = np.arange(n + 1)
    X, Y = np.meshgrid(x, y)
    
    X_tri = (X + Y) / 2
    Y_tri = (Y - X) / 2

    # 3. Visualize
    # Adjust aspect ratio to make it look like a triangle (wider than tall)
    plt.figure(figsize=(12, 4))
    
    # Use pcolormesh for the rotated grid
    # We only plot where Y >= X (the upper triangle)
    # Norm and cmap consistent with the original heatmap
    im = plt.pcolormesh(
        X_tri, Y_tri, matrix,
        cmap='RdYlBu_r',
        norm=colors.LogNorm(vmin=0.001, vmax=0.05),
        edgecolors='none'
    )
    
    # Remove the bottom half (where Y_tri < 0)
    plt.ylim(0, n / 2)
    # The X range is roughly [0, n]
    plt.xlim(0, n)
    
    plt.title(f"{region_name} (Triangular)\n{coordinates}")
    plt.colorbar(im, label="Contact Frequency", pad=0.02)
    
    # Clean up axes
    plt.gca().set_aspect('equal')
    plt.axis('off')

    # 4. Save
    output_filename = f"media/{region_name}_triangular.png"
    plt.savefig(output_filename, dpi=150, bbox_inches='tight')
    print(f"Saved to {output_filename}")
    plt.close()

def plot_compartments(file_path, region_name, coordinates):
    print(f"Processing compartments for {region_name}...")
    
    # Use 100kb for large regions to see compartments
    COMP_RESOLUTION = 100000
    uri = f"{file_path}::resolutions/{COMP_RESOLUTION}"
    
    try:
        clr = cooler.Cooler(uri)
    except Exception as e:
        print(f"Error: Could not load file. {e}")
        return

    # 1. Calculate Eigenvectors (E1) for the whole chromosome
    chrom = coordinates.split(':')[0]
    # We need the whole chromosome for eigs_cis usually, or at least a large chunk
    # bioframe.fetch_chromsizes('mm10') could be used if we had it, 
    # but we can get it from the cooler
    chromsizes = clr.chromsizes
    
    # Calculate eigenvectors for the chromosome
    print(f"  Calculating eigenvectors for {chrom} at {COMP_RESOLUTION}bp resolution...")
    cis_eigs = cooltools.eigs_cis(
        clr, 
        view_df=pd.DataFrame({'chrom': [chrom], 'start': [0], 'end': [chromsizes[chrom]], 'name': [chrom]}),
        n_eigs=1
    )
    
    # Get E1 for our specific region
    e1_track = cis_eigs[1] # [0] is the view/bins, [1] is the eigenvectors
    
    # Filter for our region
    start_bp, end_bp = map(lambda x: int(x.replace(',', '')), coordinates.split(':')[1].split('-'))
    region_bins = e1_track[(e1_track['chrom'] == chrom) & 
                           (e1_track['start'] >= start_bp) & 
                           (e1_track['end'] <= end_bp)]
    
    e1_values = region_bins['E1'].values
    
    # 2. Fetch triangular matrix for the region
    matrix = clr.matrix(balance=True).fetch(coordinates)
    n = matrix.shape[0]

    # Create coordinates for pcolormesh
    x = np.arange(n + 1)
    y = np.arange(n + 1)
    X, Y = np.meshgrid(x, y)
    X_tri = (X + Y) / 2
    Y_tri = (Y - X) / 2

    # 3. Visualize with 1D track
    fig = plt.figure(figsize=(12, 10))
    # GridSpec: 1 row for E1 (small), 1 row for triangular heatmap (large)
    gs = GridSpec(2, 1, height_ratios=[1, 3], hspace=0.05)
    
    # Top Axis: E1 Track
    ax_e1 = fig.add_subplot(gs[0])
    ax_e1.plot(np.arange(len(e1_values)), e1_values, color='black', lw=1)
    ax_e1.fill_between(np.arange(len(e1_values)), 0, e1_values, 
                       where=(e1_values > 0), color='red', alpha=0.3, interpolate=True)
    ax_e1.fill_between(np.arange(len(e1_values)), 0, e1_values, 
                       where=(e1_values < 0), color='blue', alpha=0.3, interpolate=True)
    ax_e1.axhline(0, color='grey', lw=0.5)
    ax_e1.set_ylabel("E1 (Eigenvector)")
    ax_e1.set_xlim(0, len(e1_values))
    ax_e1.set_title(f"Compartments for {region_name} ({coordinates})")
    ax_e1.set_xticks([]) # Hide x-ticks as they align with heatmap
    
    # Bottom Axis: Triangular Heatmap
    ax_heat = fig.add_subplot(gs[1])
    im = ax_heat.pcolormesh(
        X_tri, Y_tri, matrix,
        cmap='RdYlBu_r',
        norm=colors.LogNorm(vmin=0.0001, vmax=0.01), # Adjusted for 100kb
        edgecolors='none'
    )
    
    ax_heat.set_ylim(0, n / 2)
    ax_heat.set_xlim(0, n)
    ax_heat.set_aspect('equal')
    ax_heat.axis('off')
    
    # Add colorbar
    cax = fig.add_axes([0.92, 0.15, 0.02, 0.5])
    fig.colorbar(im, cax=cax, label="Contact Frequency")

    # 4. Save
    output_filename = f"media/{region_name}_triangular_track.png"
    plt.savefig(output_filename, dpi=150, bbox_inches='tight')
    print(f"Saved to {output_filename}")
    plt.close()

def plot_insulation(file_path, region_name, coordinates):
    print(f"Processing insulation for {region_name}...")
    
    uri = f"{file_path}::resolutions/{RESOLUTION}"
    try:
        clr = cooler.Cooler(uri)
    except Exception as e:
        print(f"Error: Could not load file. {e}")
        return

    matrix = clr.matrix(balance=True).fetch(coordinates)
    
    window_sizes = [25000, 50000, 100000]
    chrom = coordinates.split(':')[0]
    start_bp, end_bp = map(lambda x: int(x.replace(',', '')), coordinates.split(':')[1].split('-'))
    
    context_buffer = 500000
    fetch_start = max(0, start_bp - context_buffer)
    fetch_end = min(clr.chromsizes[chrom], end_bp + context_buffer)
    
    insulation_table = cooltools.insulation(
        clr, window_sizes, ignore_diags=2,
        view_df=pd.DataFrame({'chrom': [chrom], 'start': [fetch_start], 'end': [fetch_end], 'name': [chrom]})
    )
    
    region_ins = insulation_table[(insulation_table['start'] >= start_bp) & (insulation_table['end'] <= end_bp)]
    
    fig = plt.figure(figsize=(10, 12))
    gs = GridSpec(2, 1, height_ratios=[2, 1], hspace=0.1)
    
    ax_heat = fig.add_subplot(gs[0])
    im = ax_heat.imshow(matrix, cmap='RdYlBu_r', interpolation='none', norm=colors.LogNorm(vmin=0.001, vmax=0.05))
    ax_heat.set_title(f"Insulation & Boundaries: {region_name}")

    boundary_col = f'is_boundary_{window_sizes[1]}'
    if boundary_col in region_ins.columns:
        boundaries = region_ins[region_ins[boundary_col] == True]
        for _, row in boundaries.iterrows():
            pos = ((row['start'] + row['end']) / 2 - start_bp) / RESOLUTION
            ax_heat.axvline(pos, color='black', lw=1, ls='--')
            ax_heat.axhline(pos, color='black', lw=1, ls='--')

    ax_ins = fig.add_subplot(gs[1])
    x = np.arange(len(region_ins))
    for w in window_sizes:
        col = f'log2_insulation_score_{w}'
        if col in region_ins.columns:
            ax_ins.plot(x, region_ins[col], label=f'{w//1000}kb')
    ax_ins.set_ylabel("log2 Insulation")
    ax_ins.legend(loc='lower left', fontsize='small')
    ax_ins.axhline(0, color='grey', lw=0.5)

    output_filename = f"media/{region_name}_insulation.png"
    plt.savefig(output_filename, dpi=150, bbox_inches='tight')
    print(f"Saved to {output_filename}")
    plt.close()

def plot_dots(file_path, region_name, coordinates):
    print(f"Processing dots (loops) for {region_name}...")
    
    uri = f"{file_path}::resolutions/{RESOLUTION}"
    try:
        clr = cooler.Cooler(uri)
    except Exception as e:
        print(f"Error: Could not load file. {e}")
        return

    chrom = coordinates.split(':')[0]
    # Dots calling usually requires whole chromosome expected
    print("  Calculating expected and calling dots...")
    expected = cooltools.expected_cis(clr, view_df=pd.DataFrame({'chrom': [chrom], 'start': [0], 'end': [clr.chromsizes[chrom]], 'name': [chrom]}))
    
    # Call dots (this might be slow, usually done with GPUs or carefully)
    # For speed in this demo, we'll use a very limited set of kernels or a small region if filterable
    # But cooltools.dots usually runs on whole cooler/chrom
    dots = cooltools.dots(clr, expected=expected, nproc=1)
    
    # Filter dots for our region
    start_bp, end_bp = map(lambda x: int(x.replace(',', '')), coordinates.split(':')[1].split('-'))
    region_dots = dots[(dots['chrom1'] == chrom) & (dots['start1'] >= start_bp) & (dots['end1'] <= end_bp) &
                       (dots['chrom2'] == chrom) & (dots['start2'] >= start_bp) & (dots['end2'] <= end_bp)]
    
    matrix = clr.matrix(balance=True).fetch(coordinates)
    plt.figure(figsize=(8, 8))
    plt.imshow(matrix, cmap='RdYlBu_r', interpolation='none', norm=colors.LogNorm(vmin=0.001, vmax=0.05))
    
    # Plot dots as circles
    for _, dot in region_dots.iterrows():
        p1 = (dot['start1'] + dot['end1']) / 2 - start_bp
        p2 = (dot['start2'] + dot['end2']) / 2 - start_bp
        plt.scatter(p2 / RESOLUTION, p1 / RESOLUTION, s=40, edgecolors='black', facecolors='none', lw=1)

    plt.title(f"Dots (Loops): {region_name}")
    output_filename = f"media/{region_name}_dots.png"
    plt.savefig(output_filename, dpi=150, bbox_inches='tight')
    print(f"Saved to {output_filename}")
    plt.close()

def plot_saddle(file_path, region_name, coordinates):
    print(f"Processing saddle plot for {region_name}...")
    
    # Saddle plots usually done at coarser resolution
    SAD_RES = 100000
    uri = f"{file_path}::resolutions/{SAD_RES}"
    clr = cooler.Cooler(uri)
    chrom = coordinates.split(':')[0]
    
    # 1. E1 track
    view = pd.DataFrame({'chrom': [chrom], 'start': [0], 'end': [clr.chromsizes[chrom]], 'name': [chrom]})
    cis_eigs = cooltools.eigs_cis(clr, view_df=view, n_eigs=1)
    e1 = cis_eigs[1][['chrom', 'start', 'end', 'E1']]
    
    # 2. Expected
    expected = cooltools.expected_cis(clr, view_df=view)
    
    # 3. Saddle
    # Digitizing E1 into 50 bins
    Q_LO, Q_HI = 0.02, 0.98
    
    # Track needs to be a dataframe with chrom, start, end, and value
    # e1 already has these from eigs_cis
    
    interaction_sum, count_sum = cooltools.saddle(
        clr, expected, e1, 'cis', 50,
        qrange=(Q_LO, Q_HI),
        view_df=view
    )
    
    saddle = interaction_sum[1:-1, 1:-1] / count_sum[1:-1, 1:-1]
    plt.figure(figsize=(8, 7))
    plt.imshow(np.log10(saddle), cmap='coolwarm', vmin=-0.5, vmax=0.5)
    plt.colorbar(label="log10 (Interaction Strength / Expected)")
    plt.title(f"Saddle Plot: {chrom}\n(Ordered by E1)")
    plt.xlabel("E1 Quantile")
    plt.ylabel("E1 Quantile")
    
    output_filename = f"media/{region_name}_saddle.png"
    plt.savefig(output_filename, dpi=150, bbox_inches='tight')
    print(f"Saved to {output_filename}")
    plt.close()

if __name__ == "__main__":
    if not os.path.exists(FILE_PATH):
        print(f"STOP: I can't find the file at {FILE_PATH}")
    else:
        for name, coords in regions.items():
            if "Compartments" in name:
                plot_compartments(FILE_PATH, name, coords)
                plot_saddle(FILE_PATH, name, coords)
            else:
                plot_region(FILE_PATH, name, coords)
                plot_triangular_region(FILE_PATH, name, coords)
                plot_insulation(FILE_PATH, name, coords)
                # Note: plot_dots is computationally expensive
                # plot_dots(FILE_PATH, name, coords)