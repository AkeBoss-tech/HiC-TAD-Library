# HG-DT local reference files

Files are stored under `data/references/` (gitignored). Populate them with:

```bash
make download-references
```

Or: `bash scripts/download_hgdt_references.sh`

## What gets downloaded

| File | Source | Notes |
|------|--------|--------|
| `hg38.fa` | UCSC `hg38.fa.gz` | Primary assembly (~3 Gbp uncompressed). Resume with `curl -C -` if interrupted. |
| `gencode.v49.annotation.gtf` | GENCODE 49 | Gene/exon structure for `ReferenceContextBuilder`. |
| `GRCh38-cCREs.ELS.bed` | SCREEN Registry v4 | Candidate enhancer-like cCREs (large subset). Older “single combined” BED URLs often 404. |

## Environment overrides

- `HGDT_REF_DIR` — directory for all three (default: `<repo>/data/references`).
- `HGDT_HG38_FASTA`, `HGDT_GENCODE_GTF`, `HGDT_CCRE_BED` — individual paths.
- `HGDT_SKIP_HG38_FA=1` — skip only the hg38 FASTA step (useful when bandwidth-limited).

## Using local hg38 in the app

When `hg38.fa` is present and passes a size check, **HG-DT** reads the 1 Mb AlphaGenome windows from disk instead of the UCSC REST API. Annotations in Step 2 still use UCSC unless you extend the pipeline to query GENCODE locally.

## Citation

SCREEN data: cite Moore…Weng (*Nature* 2026) per [SCREEN downloads](https://screen.wenglab.org/downloads).
