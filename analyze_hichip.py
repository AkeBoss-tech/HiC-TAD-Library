"""HiChIP loop calling and visualisation.

Uses the real mcool file for the contact matrix.
ChIP peaks are loaded from data/raw/chip_peaks.bed if present; otherwise
synthetic CTCF-like peaks are generated at TAD boundaries (biologically
reasonable since CTCF/cohesin marks most loop anchors).

Produces:
    media/hichip_loops_sox11.png     — contact matrix + peak-anchored loops
    media/hichip_loops_mir9.png
    media/hichip_peak_profile.png    — aggregate contact profile at peaks
    media/hichip_peak_enrichment.png — per-peak enrichment scores
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import cooler

from src.hichip import call_hichip_loops, peak_contact_profile, compute_peak_enrichment


FILE_PATH = "/Users/akashdubey/Documents/CodingProjects/HiC-TAD-Library/data/raw/mouse_microc.mcool"
PEAKS_PATH = "data/raw/chip_peaks.bed"
RESOLUTION = 5_000

REGIONS = {
    'sox11_chr12': 'chr12:26,000,000-28,000,000',
    'mir9_chr13':  'chr13:83,500,000-84,500,000',
}

os.makedirs('media', exist_ok=True)


# ---------------------------------------------------------------------------
# Synthetic CTCF-like peaks at TAD boundaries
# ---------------------------------------------------------------------------

def generate_synthetic_peaks(clr, region, n_peaks=30, seed=42):
    """Generate synthetic ChIP peaks at strong insulation score minima."""
    import cooltools
    from src.tad_boundaries import parse_coordinates, make_view_df, score_boundary_prominence

    chrom, start, end = parse_coordinates(region)
    view_df = make_view_df(chrom, start, end)
    window_bp = max(5 * clr.binsize, 25_000)

    try:
        ins = cooltools.insulation(clr, [window_bp], view_df=view_df)
        ins = ins[(ins['chrom'] == chrom) & (ins['start'] >= start) & (ins['end'] <= end)]
        scored = score_boundary_prominence(ins, window_bp)
        boundaries = scored[scored['boundary_class'].isin(['strong', 'weak'])].head(n_peaks)

        if len(boundaries) == 0:
            raise ValueError("no boundaries found")

        peaks = boundaries[['chrom', 'start', 'end']].copy()
        # Widen peaks to ~300 bp centred on boundary midpoint
        mid = (peaks['start'] + peaks['end']) // 2
        peaks['start'] = (mid - 150).clip(lower=start)
        peaks['end'] = (mid + 150).clip(upper=end)
        print(f"  Generated {len(peaks)} synthetic CTCF-like peaks at TAD boundaries ({region})")
        return peaks.reset_index(drop=True)

    except Exception as exc:
        # Fallback: uniformly spaced peaks
        rng = np.random.default_rng(seed)
        n = min(n_peaks, (end - start) // clr.binsize - 2)
        starts = rng.integers(start, end - 1000, n)
        starts.sort()
        df = pd.DataFrame({
            'chrom': chrom,
            'start': starts,
            'end': starts + 300,
        })
        print(f"  Generated {len(df)} uniformly-spaced synthetic peaks ({region})")
        return df


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def plot_hichip_loops(clr, peaks_df, region, label, out_path):
    from src.tad_boundaries import parse_coordinates
    chrom, start, end = parse_coordinates(region)
    resolution = clr.binsize

    print(f"  Calling HiChIP loops for {label} …")
    try:
        loops_df = call_hichip_loops(
            clr, peaks_df, region,
            min_dist_bp=resolution * 2,
            max_dist_bp=1_500_000,
            enrichment_threshold=1.5,
            fdr_threshold=0.2,
        )
    except Exception as exc:
        print(f"  ERROR: {exc}")
        return None

    print(f"  Found {len(loops_df)} HiChIP loops in {label}")

    matrix = clr.matrix(balance=True).fetch(f"{chrom}:{start}-{end}")
    matrix_log = np.log1p(np.nan_to_num(matrix, nan=0.0))
    extent = [start / 1e6, end / 1e6, end / 1e6, start / 1e6]

    fig, ax = plt.subplots(figsize=(9, 8))
    im = ax.imshow(matrix_log, aspect='equal', origin='upper', cmap='YlOrRd', extent=extent)
    plt.colorbar(im, ax=ax, shrink=0.8, label='log(1 + balanced contact)')

    # Overlay peak positions as horizontal/vertical lines
    region_peaks = peaks_df[
        (peaks_df['chrom'] == chrom) &
        (peaks_df['start'] >= start) &
        (peaks_df['end'] <= end)
    ]
    for _, pk in region_peaks.iterrows():
        mid = (pk['start'] + pk['end']) / 2 / 1e6
        ax.axhline(mid, color='green', linewidth=0.4, alpha=0.5)
        ax.axvline(mid, color='green', linewidth=0.4, alpha=0.5)

    # Overlay loops
    for _, row in loops_df.iterrows():
        x = (row['start2'] + row['end2']) / 2 / 1e6
        y = (row['start1'] + row['end1']) / 2 / 1e6
        size = max(30, row['enrichment'] * 25)
        ax.scatter(x, y, s=size, c='#9467bd', marker='D',
                   alpha=0.8, edgecolors='white', linewidths=1, zorder=5)
        ax.scatter(y, x, s=size, c='#9467bd', marker='D',
                   alpha=0.8, edgecolors='white', linewidths=1, zorder=5)

    ax.set_xlabel('Genomic position (Mb)')
    ax.set_ylabel('Genomic position (Mb)')
    ax.set_title(
        f'HiChIP loops — {label.replace("_", " ").title()}\n'
        f'{len(region_peaks)} peaks  |  {len(loops_df)} peak-anchored loops  (5 kb)',
    )

    import matplotlib.patches as mpatches
    handles = [
        mpatches.Patch(color='green', alpha=0.5, label='ChIP peaks'),
        mpatches.Patch(color='#9467bd', label=f'HiChIP loops (n={len(loops_df)})'),
    ]
    ax.legend(handles=handles, loc='upper left', fontsize=9)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}")
    return loops_df


def plot_peak_profile(clr, peaks_df, region, label, out_path, flank_bins=25):
    print(f"  Computing peak contact profiles for {label} …")
    try:
        profiles = peak_contact_profile(clr, peaks_df, region, flank_bins=flank_bins)
    except Exception as exc:
        print(f"  ERROR: {exc}")
        return

    if profiles.shape[0] == 0:
        print("  No valid profiles.")
        return

    resolution = clr.binsize
    mean_profile = np.mean(profiles, axis=0)
    size = 2 * flank_bins + 1
    x_kb = np.arange(-flank_bins, flank_bins + 1) * resolution / 1e3

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    im = ax.imshow(profiles, aspect='auto', origin='upper', cmap='YlOrRd',
                   vmin=np.nanpercentile(profiles, 5), vmax=np.nanpercentile(profiles, 99))
    plt.colorbar(im, ax=ax, label='Balanced contact')
    ax.set_xticks(np.linspace(0, size - 1, 5))
    ax.set_xticklabels([f'{v:.0f}' for v in np.linspace(-flank_bins, flank_bins, 5) * resolution / 1e3])
    ax.set_xlabel('Distance from peak (kb)')
    ax.set_ylabel('Peak index')
    ax.set_title(f'Contact profile per peak  (n={profiles.shape[0]})')

    ax2 = axes[1]
    ax2.plot(x_kb, mean_profile, color='#9467bd', linewidth=2)
    ax2.fill_between(x_kb,
                     mean_profile - np.std(profiles, axis=0),
                     mean_profile + np.std(profiles, axis=0),
                     alpha=0.2, color='#9467bd')
    ax2.axvline(0, color='black', linestyle='--', linewidth=1)
    ax2.set_xlabel('Distance from peak centre (kb)')
    ax2.set_ylabel('Mean balanced contact')
    ax2.set_title('Mean contact profile ± 1 SD')

    plt.suptitle(f'HiChIP peak-centred profiles — {label.replace("_", " ").title()}')
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}")


def plot_peak_enrichment(clr, peaks_df, region, label, out_path):
    print(f"  Computing peak enrichment scores for {label} …")
    try:
        enr_df = compute_peak_enrichment(clr, peaks_df, region)
    except Exception as exc:
        print(f"  ERROR: {exc}")
        return

    if enr_df.empty:
        return

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))

    ax = axes[0]
    enr_sorted = enr_df.sort_values('enrichment_over_median', ascending=False)
    colors = ['#d62728' if e > 1.5 else '#1f77b4' for e in enr_sorted['enrichment_over_median']]
    ax.bar(range(len(enr_sorted)), enr_sorted['enrichment_over_median'], color=colors, edgecolor='none')
    ax.axhline(1.0, color='grey', linestyle='--', linewidth=1)
    ax.axhline(1.5, color='red', linestyle='--', linewidth=0.8, alpha=0.5, label='1.5× threshold')
    ax.set_xlabel('Peak rank')
    ax.set_ylabel('Contact enrichment over median')
    ax.set_title(f'Per-peak Hi-C contact enrichment  ({len(enr_df)} peaks)')
    ax.legend(fontsize=8)

    ax2 = axes[1]
    ax2.hist(enr_df['enrichment_over_median'], bins=20, color='#9467bd', edgecolor='white')
    ax2.axvline(1.0, color='grey', linestyle='--', linewidth=1.5, label='Median baseline')
    ax2.set_xlabel('Enrichment over median')
    ax2.set_ylabel('Number of peaks')
    ax2.set_title('Enrichment distribution')
    ax2.legend(fontsize=8)

    plt.suptitle(f'ChIP peak enrichment — {label.replace("_", " ").title()}', y=1.01)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if not os.path.exists(FILE_PATH):
        print(f"ERROR: {FILE_PATH} not found.")
        sys.exit(1)

    demo_mode = not os.path.exists(PEAKS_PATH)
    if demo_mode:
        print("NOTE: ChIP peaks not found — generating synthetic CTCF-like peaks at TAD boundaries.")
        print(f"  To use real data, place a BED file at {PEAKS_PATH}")
    else:
        print(f"Loading peaks from {PEAKS_PATH}")

    print(f"Loading {FILE_PATH} at {RESOLUTION // 1000} kb …")
    clr = cooler.Cooler(f"{FILE_PATH}::resolutions/{RESOLUTION}")

    for label, region in REGIONS.items():
        if demo_mode:
            from src.hichip import load_peaks
            peaks_df = generate_synthetic_peaks(clr, region)
        else:
            from src.hichip import load_peaks
            from src.tad_boundaries import parse_coordinates
            chrom, _, _ = parse_coordinates(region)
            peaks_df = load_peaks(PEAKS_PATH, chrom=chrom)

        plot_hichip_loops(clr, peaks_df, region, label, f"media/hichip_loops_{label}.png")
        plot_peak_profile(clr, peaks_df, region, label, f"media/hichip_peak_profile_{label}.png")
        plot_peak_enrichment(clr, peaks_df, region, label, f"media/hichip_peak_enrichment_{label}.png")

    print("\nDone. Outputs in media/")
    print("\nResult interpretation:")
    print("  HiChIP loops  : Purple diamonds mark contact pixels where BOTH anchors")
    print("    overlap a ChIP peak.  Green lines show peak positions. In real HiChIP")
    print("    (e.g. SMC1 or H3K27ac), peaks are enriched at loop anchors and enhancers.")
    print("  Peak profiles : Each row = contact profile for one peak.  A peak with")
    print("    strong looping shows elevated contact at a specific offset (the loop partner).")
    print("  Enrichment    : Peaks enriched >1.5× over background are likely at loop")
    print("    anchors or active regulatory elements.")


if __name__ == '__main__':
    main()
