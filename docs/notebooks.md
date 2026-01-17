# Notebooks Documentation

This directory contains Jupyter notebooks and scripts for analyzing Hi-C data.

## Overview

The analysis workflow is divided into three main stages:
1. **Data Fetching**: `00_fetch_data.py`
2. **Exploration**: `01_explore.ipynb`
3. **TAD Analysis**: `02_tad_analysis.ipynb`
4. **Compartment Analysis**: `03_compartments.ipynb`

---

## Detailed Descriptions

### [00_fetch_data.py](../notebooks/00_fetch_data.py)
**Purpose**: download test data for the library.
- Fetches `test.mcool` (HFF Micro-C data) from a remote source if not present.
- Saves the file to `data/raw/`.

**Usage**:
```bash
python notebooks/00_fetch_data.py
```

### [01_explore.ipynb](../notebooks/01_explore.ipynb)
**Purpose**: Load and visualize Contact Matrices.
- **Key Features**:
    - Demonstrates how to use `src.loaders.get_cooler` to load `.mcool` files.
    - Visualizes the contact matrix for a specific chromosome (e.g., `chr17`) using `matplotlib`.
    - Applies `log1p` transformation for better visibility of contact frequencies.

### [02_tad_analysis.ipynb](../notebooks/02_tad_analysis.ipynb)
**Purpose**: Identify and visualize Topologically Associating Domains (TADs).
- **Key Features**:
    - Calculates the insulation score using `cooltools.insulation`.
    - Calls TAD boundaries based on insulation minima.
    - Visualizes the contact matrix with overlaid TAD boundaries.

### [03_compartments.ipynb](../notebooks/03_compartments.ipynb)
**Purpose**: A/B Compartment Analysis.
- **Key Features**:
    - Performs eigenvector decomposition using `cooltools.eigs_cis` to separate active (A) and inactive (B) compartments.
    - Visualizes the first eigenvector (PC1) alongside the contact map.
    - *Note*: Requires `test.mcool` in `data/raw/` (fixed file path issue).
