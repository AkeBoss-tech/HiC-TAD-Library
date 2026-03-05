"""SNV and small indel analysis from BAM files.

If data/raw/aligned.bam and data/raw/mm10.fa are present, runs the real
pileup-based caller.  Otherwise, generates a synthetic SNV dataset that
demonstrates the analysis pipeline (VAF spectrum, Ti/Tv ratio, etc.).

Produces:
    media/variants_vaf_histogram.png
    media/variants_spectrum.png
    media/variants_ti_tv.png
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

from src.variants import annotate_ti_tv

BAM_PATH   = "data/raw/aligned.bam"
REF_FASTA  = "data/raw/mm10.fa"
REGION     = "chr12:26,000,000-28,000,000"

os.makedirs('media', exist_ok=True)


# ---------------------------------------------------------------------------
# Synthetic demo data (used when BAM is absent)
# ---------------------------------------------------------------------------

def generate_synthetic_variants(region=REGION, n_snv=800, n_indel=80, seed=42):
    """Generate realistic synthetic SNV/indel calls for demonstration."""
    rng = np.random.default_rng(seed)
    from src.tad_boundaries import parse_coordinates
    chrom, start, end = parse_coordinates(region)

    # SNVs
    pos = rng.integers(start, end, n_snv)
    bases = list('ACGT')
    refs = rng.choice(bases, n_snv)
    alts = np.array([rng.choice([b for b in bases if b != r]) for r in refs])

    # Bimodal VAF: germline het (~0.5) and somatic mosaic (~0.1–0.3)
    vaf_germline = rng.normal(0.5, 0.04, n_snv // 2)
    vaf_somatic = rng.beta(2, 8, n_snv - n_snv // 2) * 0.4 + 0.05
    vafs = np.clip(np.concatenate([vaf_germline, vaf_somatic]), 0.05, 0.95)
    rng.shuffle(vafs)

    depths = rng.integers(15, 120, n_snv)
    alt_counts = (vafs * depths).astype(int)
    quals = (depths * vafs * 0.9).astype(int)

    snv_df = pd.DataFrame({
        'chrom': chrom,
        'pos': pos,
        'ref': refs,
        'alt': alts,
        'depth': depths,
        'alt_count': alt_counts,
        'vaf': vafs,
        'strand_bias': rng.uniform(0, 0.3, n_snv),
        'strand_bias_flag': rng.choice([True, False], n_snv, p=[0.05, 0.95]),
        'qual': quals,
        'variant_type': 'snv',
    })

    # Indels
    pos_i = rng.integers(start, end, n_indel)
    indel_df = pd.DataFrame({
        'chrom': chrom,
        'pos': pos_i,
        'ref': ['A'] * n_indel,
        'alt': [rng.choice(['AG', 'AT', '-', 'ATCG']) for _ in range(n_indel)],
        'depth': rng.integers(15, 100, n_indel),
        'alt_count': rng.integers(5, 50, n_indel),
        'vaf': rng.uniform(0.1, 0.6, n_indel),
        'strand_bias': rng.uniform(0, 0.3, n_indel),
        'strand_bias_flag': rng.choice([True, False], n_indel, p=[0.08, 0.92]),
        'qual': rng.integers(10, 80, n_indel),
        'variant_type': 'indel',
    })

    return pd.concat([snv_df, indel_df], ignore_index=True).sort_values('pos').reset_index(drop=True)


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def plot_vaf_histogram(df, out_path):
    snv = df[df['variant_type'] == 'snv']
    indel = df[df['variant_type'] == 'indel']

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    ax.hist(snv['vaf'].dropna(), bins=40, color='#1f77b4', edgecolor='white', alpha=0.8)
    ax.axvline(0.5, color='red', linestyle='--', linewidth=1.2, label='Germline het (0.5)')
    ax.set_xlabel('Variant allele frequency (VAF)')
    ax.set_ylabel('Number of SNVs')
    ax.set_title(f'SNV VAF distribution  (n={len(snv):,})')
    ax.legend(fontsize=9)

    ax2 = axes[1]
    ax2.hist(indel['vaf'].dropna(), bins=20, color='#ff7f0e', edgecolor='white', alpha=0.8)
    ax2.set_xlabel('Variant allele frequency (VAF)')
    ax2.set_ylabel('Number of indels')
    ax2.set_title(f'Indel VAF distribution  (n={len(indel):,})')

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}")


def plot_mutation_spectrum(df, out_path):
    """Plot the 6-class substitution spectrum (SBS context)."""
    snv = df[(df['variant_type'] == 'snv')].copy()
    snv = annotate_ti_tv(snv)

    # 6 substitution classes (on pyrimidine strand)
    def normalise_sub(ref, alt):
        purines = {'A', 'G'}
        if ref in purines:
            comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
            return comp.get(ref, ref), comp.get(alt, alt)
        return ref, alt

    sub_counts = {}
    for _, row in snv.iterrows():
        r, a = str(row['ref']).upper(), str(row['alt']).upper()
        if len(r) == 1 and len(a) == 1:
            r2, a2 = normalise_sub(r, a)
            key = f'{r2}>{a2}'
            sub_counts[key] = sub_counts.get(key, 0) + 1

    expected_keys = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
    colors = ['#1f77b4', '#000000', '#d62728', '#aec7e8', '#2ca02c', '#ff7f0e']

    counts = [sub_counts.get(k, 0) for k in expected_keys]
    total = sum(counts) or 1
    freqs = [c / total for c in counts]

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(expected_keys, freqs, color=colors, edgecolor='white')
    ax.set_ylabel('Fraction of SNVs')
    ax.set_title('Substitution spectrum (SBS)\nC>T transitions dominate in normal cells')
    ax.set_ylim(0, max(freqs) * 1.3)

    for bar, count in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.005,
                f'{count:,}', ha='center', va='bottom', fontsize=9)

    # Ti/Tv ratio
    ti_count = snv['ti_tv'].eq('transition').sum() if 'ti_tv' in snv.columns else 0
    tv_count = snv['ti_tv'].eq('transversion').sum() if 'ti_tv' in snv.columns else 0
    ti_tv_ratio = ti_count / max(tv_count, 1)
    ax.text(0.98, 0.95, f'Ti/Tv ratio: {ti_tv_ratio:.2f}',
            transform=ax.transAxes, ha='right', va='top', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}")


def plot_ti_tv_summary(df, out_path):
    snv = annotate_ti_tv(df[df['variant_type'] == 'snv'])

    ti = (snv['ti_tv'] == 'transition').sum()
    tv = (snv['ti_tv'] == 'transversion').sum()
    ratio = ti / max(tv, 1)

    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    ax = axes[0]
    ax.pie([ti, tv], labels=[f'Transitions\n(n={ti:,})', f'Transversions\n(n={tv:,})'],
           colors=['#1f77b4', '#d62728'], autopct='%1.1f%%', startangle=90)
    ax.set_title('Transition vs Transversion')

    ax2 = axes[1]
    ax2.barh(['Ti/Tv ratio'], [ratio], color='#2ca02c', height=0.4)
    ax2.axvline(2.0, color='orange', linestyle='--', label='Expected germline (~2.0)')
    ax2.axvline(3.3, color='red', linestyle='--', label='Expected exome (~3.3)')
    ax2.set_xlabel('Ti/Tv ratio')
    ax2.set_title('Ti/Tv quality metric\n(low ratio suggests excess noise)')
    ax2.set_xlim(0, max(4, ratio + 0.5))
    ax2.legend(fontsize=8)
    ax2.text(ratio + 0.05, 0, f'{ratio:.2f}', va='center', fontsize=10)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    demo_mode = not (os.path.exists(BAM_PATH) and os.path.exists(REF_FASTA))

    if demo_mode:
        print("NOTE: BAM/FASTA not found — generating synthetic variant data for demonstration.")
        print(f"  To use real data, place aligned BAM at {BAM_PATH} and FASTA at {REF_FASTA}")
        df = generate_synthetic_variants()
        print(f"  Generated {len(df)} synthetic variants ({REGION})")
    else:
        from src.variants import call_snvs
        print(f"Calling variants in {REGION} from {BAM_PATH} …")
        try:
            df = call_snvs(BAM_PATH, REGION, ref_fasta=REF_FASTA)
            print(f"  Called {len(df)} variants")
        except Exception as exc:
            print(f"  ERROR: {exc}")
            sys.exit(1)

    if df.empty:
        print("No variants to plot. Exiting.")
        sys.exit(0)

    df = annotate_ti_tv(df)

    plot_vaf_histogram(df, "media/variants_vaf_histogram.png")
    plot_mutation_spectrum(df, "media/variants_spectrum.png")
    plot_ti_tv_summary(df, "media/variants_ti_tv.png")

    # Summary
    snv = df[df['variant_type'] == 'snv']
    indel = df[df['variant_type'] == 'indel']
    ti = (snv['ti_tv'] == 'transition').sum() if 'ti_tv' in snv.columns else 0
    tv = (snv['ti_tv'] == 'transversion').sum() if 'ti_tv' in snv.columns else 0

    print("\nDone. Outputs in media/")
    print("\nVariant summary:")
    print(f"  SNVs:   {len(snv):>6,}")
    print(f"  Indels: {len(indel):>6,}")
    print(f"  Ti:     {ti:>6,}  |  Tv: {tv:,}  |  Ti/Tv: {ti/max(tv,1):.2f}")
    if not demo_mode:
        print(f"  Median depth:   {df['depth'].median():.0f}×")
        print(f"  Median VAF:     {df['vaf'].median():.2f}")
    print("\nResult interpretation:")
    print("  VAF histogram : Germline het variants cluster at VAF~0.5; somatic mosaic")
    print("    variants appear at low VAF (0.05–0.3). Two-peak distribution is expected")
    print("    in normal tissue.")
    print("  Mutation spectrum : C>T transitions dominate in normal genomes (deamination).")
    print("    APOBEC-driven tumours show C>T/G in TC context; UV damage shows CC>TT.")
    print("  Ti/Tv ratio : ~2.0 expected genome-wide; exomes ~3.3. Values <1.5 suggest")
    print("    sequencing noise or relaxed filters.")


if __name__ == '__main__':
    main()
