# Multi-Modal Deletion Analysis with AlphaGenome

This document describes the enhanced deletion analysis pipeline that now predicts **contact maps**, **TF binding (ChIP-seq)**, and **ATAC-seq accessibility** for mouse genomic deletions.

## What's New

### Additional Predictions
Previously, the deletion analysis only predicted contact maps. Now it predicts:

1. **Contact Maps** - 3D chromatin architecture (unchanged)
2. **TF Binding (CHIP_TF)** - Transcription factor binding predictions across multiple TFs
3. **ATAC-seq (ATAC)** - Chromatin accessibility predictions showing open regulatory regions

### New Visualizations

#### 1. Multi-Modal Overview (`*_multimodal.png`)
A comprehensive figure combining all three data types:
- Row 1: Contact maps (WT | Deletion | Difference)
- Row 2: Mean TF binding tracks (WT vs Deletion)
- Row 3: Mean ATAC-seq tracks (WT vs Deletion)

This gives you an integrated view of how the deletion affects:
- 3D chromatin structure
- Transcription factor binding
- Chromatin accessibility

#### 2. Detailed 1D Tracks (`*_1d_tracks.png`)
Four-panel figure showing:
1. **TF Binding Tracks** - Multiple TF predictions (up to 5 TFs shown)
   - Solid lines = Wild-type
   - Dashed lines = After deletion
2. **TF Binding Difference** - Mean change (red = gained, blue = lost)
3. **ATAC-seq Tracks** - Multiple accessibility predictions
   - Solid lines = Wild-type
   - Dashed lines = After deletion
4. **ATAC Difference** - Mean change (red = gained, blue = lost)

### Interpretation Guide

#### TF Binding Changes
- **Gained binding (red)**: Deletion may expose new TF binding sites or create favorable chromatin context
- **Lost binding (blue)**: Deletion removes TF binding sites or disrupts chromatin accessibility
- **Look for**: Changes at the deletion boundaries that may indicate altered regulatory element positioning

#### ATAC-seq Changes
- **Gained accessibility (red)**: Deletion may open up chromatin or bring enhancers closer
- **Lost accessibility (blue)**: Deletion removes accessible regions or disrupts chromatin loops
- **Look for**: Changes in peak positions that correlate with TAD boundary shifts

#### Combined Interpretation
Compare across modalities to understand mechanism:
- Contact map changes + TF binding loss → Deletion disrupts loop anchoring
- Contact map changes + ATAC changes → Deletion affects enhancer-promoter contacts
- All three changing → Major chromatin reorganization

## Running the Analysis

```bash
# Run the enhanced deletion analysis
python run_mouse_deletion.py

# This will generate:
# - All original contact map visualizations
# - New TF binding and ATAC-seq tracks
# - Multi-modal overview figures
# - Updated HTML report with all modalities
```

## Output Files

### Per Cell Type (for each deletion region):
- `mouse_deletion_*_multimodal.png` - Multi-modal overview
- `mouse_deletion_*_1d_tracks.png` - Detailed TF & ATAC tracks
- `mouse_deletion_*_triangle.png` - Triangle TAD view (contact maps only)
- `mouse_deletion_*.png` - Square contact map with gene track
- `mouse_deletion_*_extra.png` - Log2 ratio, virtual 4C, P(s) curves

### Combined Reports:
- `analysis_report.html` - Interactive HTML report with all visualizations

## Technical Details

### AlphaGenome Output Types Used
```python
from alphagenome.models import dna_client

requested_outputs = {
    dna_client.OutputType.CONTACT_MAPS,  # 3D contacts
    dna_client.OutputType.CHIP_TF,       # TF binding
    dna_client.OutputType.ATAC,          # Accessibility
}
```

### Data Shapes
- **Contact maps**: `(n_bins, n_bins, n_cell_types)` - 2D symmetric matrix
- **TF binding**: `(n_bins, n_TFs)` - Multiple TF tracks per position
- **ATAC-seq**: `(n_bins, n_tracks)` - Multiple accessibility tracks

### Resolution
- Genomic window: 1 MB (2^20 bp) centered on deletion (AlphaGenome maximum)
- Contact map resolution: Model-dependent (typically ~2048 bp)
- Track resolution: Model-dependent (typically ~128-512 bp)

## Example Use Cases

### 1. TAD Boundary Deletion
**Expected effects:**
- Contact maps: TAD fusion, loss of insulation
- TF binding: Loss of CTCF at boundary
- ATAC: Reduced accessibility at boundary region

### 2. Enhancer Deletion
**Expected effects:**
- Contact maps: Loss of specific loop
- TF binding: Loss of enhancer-binding TFs
- ATAC: Major loss of accessibility at deletion site

### 3. Regulatory Element Deletion
**Expected effects:**
- Contact maps: Subtle changes in local contacts
- TF binding: Specific TF peaks lost
- ATAC: Localized accessibility changes

## Limitations

1. **Cell type availability**: Only 8 mouse contact map tracks in AlphaGenome training data
   - No inner-ear specific cell types available
   - Using olfactory receptor cells (CL:0000207) as proxy
   - Using embryonic stem cells (EFO:0004038) as reference

2. **Window size**: Limited to 1 MB (AlphaGenome API constraint)
   - Full 2-TAD windows may exceed this limit
   - Analysis focuses on region immediately around deletion

3. **Track interpretation**:
   - AlphaGenome predicts multiple TFs but doesn't identify which specific TFs
   - Track indices (TF1, TF2, etc.) are ordered by the model
   - ATAC tracks represent different cell types or conditions

## Next Steps

To add more output types (see AlphaGenome architecture):
- `RNA_SEQ` - Gene expression changes
- `CHIP_HISTONE` - Histone modification predictions
- `DNASE` - DNase I hypersensitivity
- `CAGE` / `PROCAP` - Transcription start sites
- `SPLICE_JUNCTIONS` - Splicing changes

See `externals/alphagenome/src/alphagenome/models/dna_output.py` for all available output types.

## References

- AlphaGenome paper: https://www.nature.com/articles/s41586-023-06484-x
- AlphaGenome SDK: https://github.com/deepmind/alphagenome
- 4DNucleome data: https://data.4dnucleome.org/
