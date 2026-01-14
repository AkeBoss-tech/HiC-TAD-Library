# Hi-C 3D Genome Organization Research

This repository contains code and analysis pipelines for investigating 3D Genome Organization, focusing on Topologically Associating Domains (TADs) and chromosomal territories.

## Project Structure

```text
.
├── data/
│   ├── raw/                # Original .cool and .mcool files
│   ├── processed/          # Derived data (insulation scores, etc.)
│   └── external/           # Reference genomes and other external resources
├── notebooks/              # Jupyter notebooks for interactive analysis
├── src/                    # Reusable Python source code
│   ├── loaders.py          # Data loading utilities
│   └── analysis.py         # Analysis functions
├── environment.yml         # Conda environment definition
└── README.md
```

## Getting Started

1.  **Create the environment:**
    ```bash
    conda env create -f environment.yml
    conda activate hic-analysis
    ```

2.  **Download Data:**
    - Place `.mcool` files in `data/raw/`.

3.  **Run Notebooks:**
    - Start with `notebooks/01_explore.ipynb` to familiarize yourself with the data structures.
