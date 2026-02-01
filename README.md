# HiC-TAD-Library

A Python library for visualizing and analyzing Hi-C data, focusing on Topologically Associating Domains (TADs), compartments, and insulation scores.

## Installation

You can install the required dependencies using the provided `Makefile`:

```bash
make install
```

Required packages include: `cooler`, `cooltools`, `bioframe`, `matplotlib`, `pandas`, and `numpy`.

## Usage

Generated visualizations are saved in the `media/` directory.

### Quick Start
To install dependencies and run all visualizations:
```bash
make setup
```

### Individual Commands
- **Run Visualizations**: `make run`
- **Clean Output**: `make clean`

## Visualizations

### 1. High-Resolution Heatmaps
Standard square heatmaps for identifying TADs and local features.

![Sox11 Heatmap](media/Sox11_Chr12_heatmap.png)

### 2. Triangular Heatmaps
Rotated 45-degree views, standard for visualizing domains.

![Sox11 Triangular](media/Sox11_Chr12_triangular.png)

### 3. Insulation Scores & Boundaries
Detecting structural boundaries using sliding window insulation scores.

![Sox11 Insulation](media/Sox11_Chr12_insulation.png)

### 4. Compartments (E1 Track)
Large-scale genomic organization (A/B compartments) visualized with eigenvectors.

![Compartments E1](media/Compartments_Chr2_triangular_track.png)

### 5. Saddle Plots
Analyzing compartment-dependent interaction preferences.

![Saddle Plot](media/Compartments_Chr2_saddle.png)

---
*Data sourced from [4DNucleome](https://data.4dnucleome.org/)*
