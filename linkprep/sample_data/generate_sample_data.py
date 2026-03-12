"""
Generate synthetic LinkPrep sample data.

Produces four files in sample_data/output/:
  linkprep_sample.cool          — contact matrix with TAD structure
  linkprep_sample.pairs.gz      — 4DN pairs file (ligation events)
  linkprep_sample.vcf           — heterozygous SNPs for phasing
  linkprep_sample_coverage.bed  — per-bin read-depth

Run:
    python sample_data/generate_sample_data.py
"""

import gzip
import os
import sys

import numpy as np
import pandas as pd

# ensure we can import cooler even when run from inside sample_data/
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

try:
    import cooler
except ImportError:
    sys.exit("cooler is required.  Install with: pip install cooler")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

OUT_DIR = os.path.join(os.path.dirname(__file__), "output")
os.makedirs(OUT_DIR, exist_ok=True)

COOL_PATH    = os.path.join(OUT_DIR, "linkprep_sample.cool")
PAIRS_PATH   = os.path.join(OUT_DIR, "linkprep_sample.pairs.gz")
VCF_PATH     = os.path.join(OUT_DIR, "linkprep_sample.vcf")
COV_PATH     = os.path.join(OUT_DIR, "linkprep_sample_coverage.bed")

# ---------------------------------------------------------------------------
# Genome layout (tiny synthetic chromosomes for speed)
# ---------------------------------------------------------------------------

CHROM_SIZES = {
    "chr12": 120_000_000,
    "chr13": 100_000_000,
}
BIN_SIZE = 10_000   # 10 kb

RNG = np.random.default_rng(42)

# ---------------------------------------------------------------------------
# Helper: TAD-structured contact matrix
# ---------------------------------------------------------------------------

def _tad_block_matrix(n_bins: int, tad_size: int = 20) -> np.ndarray:
    """
    Build an n_bins × n_bins symmetric contact matrix with block-TAD structure.

    Contacts decay with genomic distance; bins inside the same TAD get an
    extra 3× enrichment.  A few bright pixels simulate chromatin loops.
    """
    mat = np.zeros((n_bins, n_bins), dtype=np.float32)

    # Distance-decay background
    for d in range(n_bins):
        val = np.exp(-d / (tad_size * 2.0)) * 100 + RNG.poisson(2, n_bins - d)
        val = val.astype(np.float32)
        idx = np.arange(n_bins - d)
        mat[idx, idx + d] += val
        if d > 0:
            mat[idx + d, idx] += val

    # TAD enrichment blocks
    n_tads = n_bins // tad_size
    for t in range(n_tads):
        s, e = t * tad_size, min((t + 1) * tad_size, n_bins)
        mat[s:e, s:e] *= 3.0

    # Synthetic loops: a few bright off-diagonal pixels
    n_loops = max(1, n_bins // (tad_size * 2))
    for _ in range(n_loops):
        i = RNG.integers(0, n_bins - tad_size)
        j = i + tad_size + RNG.integers(0, tad_size)
        j = min(j, n_bins - 1)
        boost = RNG.integers(50, 150)
        mat[i, j] += boost
        mat[j, i] += boost

    np.fill_diagonal(mat, 0)
    return mat


# ---------------------------------------------------------------------------
# 1. Contact matrix (.cool)
# ---------------------------------------------------------------------------

def make_cool():
    print("Generating contact matrix (.cool) …")

    bins_list = []
    for chrom, length in CHROM_SIZES.items():
        starts = np.arange(0, length, BIN_SIZE)
        ends   = np.minimum(starts + BIN_SIZE, length)
        df = pd.DataFrame({"chrom": chrom, "start": starts, "end": ends})
        bins_list.append(df)
    bins_df = pd.concat(bins_list, ignore_index=True)

    pixels_list = []
    offset = 0
    for chrom, length in CHROM_SIZES.items():
        n = int(np.ceil(length / BIN_SIZE))
        mat = _tad_block_matrix(n, tad_size=20)

        rows, cols = np.triu_indices(n, k=1)
        counts = mat[rows, cols].astype(int)
        mask = counts > 0

        pixels_list.append(pd.DataFrame({
            "bin1_id": rows[mask] + offset,
            "bin2_id": cols[mask] + offset,
            "count":   counts[mask],
        }))
        offset += n

    pixels_df = pd.concat(pixels_list, ignore_index=True)

    cooler.create_cooler(
        COOL_PATH,
        bins=bins_df,
        pixels=pixels_df,
        dtypes={"count": np.int32},
        ordered=True,
    )

    print(f"  Saved → {COOL_PATH}")
    print("  (Raw counts — no ICE balancing applied to sample data.)")
    print("  For real data, run: cooler balance <file.cool>")


# ---------------------------------------------------------------------------
# 2. Pairs file (.pairs.gz)
# ---------------------------------------------------------------------------

PAIRS_HEADER = """\
## pairs format v1.0
#sorted: chr1-chr2-pos1-pos2
#shape: upper triangle
#chromosomes: {chroms}
#genome_assembly: mm10
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
"""

def make_pairs():
    print("Generating pairs file (.pairs.gz) …")

    chroms = list(CHROM_SIZES.keys())
    n_pairs = 500_000

    # Mostly cis (same chromosome), short-range
    cis_frac = 0.85
    n_cis   = int(n_pairs * cis_frac)
    n_trans = n_pairs - n_cis

    records = []

    # Cis pairs — distance-decay distribution
    for chrom, length in CHROM_SIZES.items():
        n = n_cis // len(CHROM_SIZES)
        pos1 = RNG.integers(1000, length - 1000, size=n)
        # distance drawn from log-normal to mimic proximity ligation
        dist  = np.exp(RNG.normal(np.log(50_000), 1.5, size=n)).astype(int)
        dist  = np.clip(dist, 500, length // 2)
        pos2  = np.clip(pos1 + dist, 1, length - 1)
        strand1 = RNG.choice(["+", "-"], size=n)
        strand2 = RNG.choice(["+", "-"], size=n)
        for i in range(n):
            p1, p2 = sorted([pos1[i], pos2[i]])
            records.append((f"read{len(records)}", chrom, p1, chrom, p2,
                            strand1[i], strand2[i], "UU"))

    # Trans (inter-chromosomal) pairs — simulate translocations + noise
    for _ in range(n_trans):
        c1, c2 = RNG.choice(chroms, size=2, replace=False)
        p1 = RNG.integers(1000, CHROM_SIZES[c1] - 1000)
        p2 = RNG.integers(1000, CHROM_SIZES[c2] - 1000)
        s1 = RNG.choice(["+", "-"])
        s2 = RNG.choice(["+", "-"])
        records.append((f"read{len(records)}", c1, p1, c2, p2, s1, s2, "UU"))

    # Inject a cluster of trans pairs simulating a translocation hotspot
    sv_pos_c1 = 26_500_000   # chr12 breakpoint
    sv_pos_c2 = 84_000_000   # chr13 breakpoint
    for _ in range(2_000):
        p1 = int(RNG.normal(sv_pos_c1, 50_000))
        p2 = int(RNG.normal(sv_pos_c2, 50_000))
        records.append((f"read{len(records)}", "chr12", p1, "chr13", p2,
                        "+", "-", "UU"))

    chrom_order = {c: i for i, c in enumerate(chroms)}
    records.sort(key=lambda r: (chrom_order.get(r[1], 99),
                                chrom_order.get(r[3], 99), r[2], r[4]))

    header = PAIRS_HEADER.format(chroms=" ".join(chroms))
    with gzip.open(PAIRS_PATH, "wt") as fh:
        fh.write(header)
        for rec in records:
            fh.write("\t".join(str(x) for x in rec) + "\n")

    print(f"  Saved → {PAIRS_PATH}  ({len(records):,} pairs)")


# ---------------------------------------------------------------------------
# 3. VCF (heterozygous SNPs for phasing)
# ---------------------------------------------------------------------------

VCF_HEADER = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depth">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
"""

BASES = list("ACGT")

def make_vcf():
    print("Generating VCF (SNPs for phasing) …")

    records = []
    snp_id  = 0
    for chrom, length in CHROM_SIZES.items():
        # ~1 SNP per 10 kb on average
        positions = sorted(RNG.integers(1000, length - 1000,
                                        size=length // 10_000))
        for pos in positions:
            ref = RNG.choice(BASES)
            alt = RNG.choice([b for b in BASES if b != ref])
            dp  = RNG.integers(20, 60)
            ad0 = RNG.integers(dp // 3, dp * 2 // 3)
            ad1 = dp - ad0
            records.append(
                f"{chrom}\t{pos}\trs{snp_id}\t{ref}\t{alt}\t"
                f"50\tPASS\t.\tGT:AD:DP\t0/1:{ad0},{ad1}:{dp}"
            )
            snp_id += 1

    with open(VCF_PATH, "w") as fh:
        fh.write(VCF_HEADER)
        fh.write("\n".join(records) + "\n")

    print(f"  Saved → {VCF_PATH}  ({len(records):,} SNPs)")


# ---------------------------------------------------------------------------
# 4. Coverage BED
# ---------------------------------------------------------------------------

def make_coverage():
    print("Generating coverage BED …")

    rows = []
    for chrom, length in CHROM_SIZES.items():
        n = int(np.ceil(length / BIN_SIZE))
        starts = np.arange(n) * BIN_SIZE
        ends   = np.minimum(starts + BIN_SIZE, length)
        cov    = RNG.poisson(30, size=n).astype(float)

        # Simulate a heterozygous deletion (~0.5× coverage) on chr12 26–27 Mb
        if chrom == "chr12":
            del_s = 26_000_000 // BIN_SIZE
            del_e = 27_000_000 // BIN_SIZE
            cov[del_s:del_e] *= 0.5

        # Simulate a duplication (~2× coverage) on chr13 83–84 Mb
        if chrom == "chr13":
            dup_s = 83_000_000 // BIN_SIZE
            dup_e = 84_000_000 // BIN_SIZE
            cov[dup_s:dup_e] *= 2.0

        for s, e, c in zip(starts, ends, cov):
            rows.append((chrom, s, e, round(c, 2)))

    df = pd.DataFrame(rows, columns=["chrom", "start", "end", "coverage"])
    df.to_csv(COV_PATH, sep="\t", index=False, header=False)
    print(f"  Saved → {COV_PATH}  ({len(df):,} bins)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    make_cool()
    make_pairs()
    make_vcf()
    make_coverage()
    print("\nAll sample data generated successfully.")
    print(f"Output directory: {OUT_DIR}")
