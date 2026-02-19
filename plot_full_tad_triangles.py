#!/usr/bin/env python3
"""
Plot the full 2-TAD regions for Ed and Jingyun as triangle contact maps
using real Mouse Micro-C data (mm10, 5 kb resolution).

Ed's TAD pair:     chr12:26,440,002–28,560,000  (2.12 Mb, 424 bins @ 5 kb)
Jingyun's TAD pair: chr13:81,760,002–85,200,000  (3.44 Mb, 688 bins @ 5 kb)

The insulator deletion sites are marked on each plot.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cooler

FILE     = 'data/raw/mouse_microc.mcool'
RES      = 5000   # 5 kb

# ── Region definitions ────────────────────────────────────────────────────────
REGIONS = [
    dict(
        name        = "Edward — chr12 full TAD pair",
        query       = "chr12:26440002-28560000",
        chrom       = "chr12",
        start       = 26_440_002,
        end         = 28_560_000,
        ins_start   = 27_333_532,
        ins_end     = 27_336_455,
        ins_label   = "Insulator\n(27.33 Mb)",
        safe        = "Edward_chr12_full_tad",
    ),
    dict(
        name        = "Jingyun — chr13 full TAD pair",
        query       = "chr13:81760002-85200000",
        chrom       = "chr13",
        start       = 81_760_002,
        end         = 85_200_000,
        ins_start   = 83_739_797,
        ins_end     = 83_745_138,
        ins_label   = "Insulator\n(83.74 Mb)",
        safe        = "Jingyun_chr13_full_tad",
    ),
]

# ── Load cooler ───────────────────────────────────────────────────────────────
print(f"Loading {FILE} @ {RES} bp ...")
clr = cooler.Cooler(f'{FILE}::resolutions/{RES}')


# ── Triangle rotation ─────────────────────────────────────────────────────────
def rotate45(matrix):
    """
    Vectorised 45° rotation for classic triangle Hi-C plot.

    element matrix[i,j] → output[j-i, i+j]
      x axis (columns, 0..2n): genomic midpoint of contact
      y axis (rows,    0..n):  genomic distance / 2  (0 = diagonal, n = max)

    Returns (n, 2n) float array; NaN where no data.
    """
    n = matrix.shape[0]
    out = np.full((n, 2 * n), np.nan)
    rows, cols = np.triu_indices(n)
    out[cols - rows, rows + cols] = matrix[rows, cols]
    return out


# ── Per-region plot ───────────────────────────────────────────────────────────
def plot_triangle(region, matrix):
    n   = matrix.shape[0]
    res = RES
    start, end = region['start'], region['end']
    span_mb    = (end - start) / 1e6
    y_max      = (end - start) / 2   # in bp

    tri = rotate45(matrix)

    # Replace NaN with 0 for display but keep masked for colour
    masked = np.ma.masked_invalid(tri)

    # ── colour scale ─────────────────────────────────────────────────────────
    # Use log-norm; clip low tail of non-zero values
    vals = matrix[~np.isnan(matrix) & (matrix > 0)]
    vmin = float(np.percentile(vals, 5))   if len(vals) else 1e-4
    vmax = float(np.percentile(vals, 99))  if len(vals) else 1.0

    cmap = plt.get_cmap('YlOrRd').copy()
    cmap.set_bad('white')
    norm = mcolors.LogNorm(vmin=max(vmin, 1e-5), vmax=vmax)

    # ── figure ────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(16, 6))
    fig.patch.set_facecolor('#fafafa')
    ax.set_facecolor('white')

    extent = [start, end, 0, y_max]
    im = ax.imshow(
        masked,
        cmap=cmap, norm=norm,
        aspect='auto',
        extent=extent,
        origin='lower',
        interpolation='nearest',
    )

    # ── insulator line ────────────────────────────────────────────────────────
    ins_mid = (region['ins_start'] + region['ins_end']) / 2
    ax.axvline(ins_mid, color='#00bfff', lw=2.5, ls='--',
               label=f"Insulator ({ins_mid/1e6:.3f} Mb)")
    ax.text(ins_mid, y_max * 0.96, region['ins_label'],
            ha='center', va='top', fontsize=8, color='#007acc',
            fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='#007acc',
                      alpha=0.85, lw=1))

    # ── TAD boundary guide lines (light grey every 0.5 Mb) ───────────────────
    mb_tick = 500_000
    for pos in range(int(start / mb_tick) * mb_tick, end + mb_tick, mb_tick):
        if start < pos < end:
            ax.axvline(pos, color='#ccc', lw=0.5, ls=':')

    # ── axes & labels ─────────────────────────────────────────────────────────
    ax.set_xlim(start, end)
    ax.set_ylim(0, y_max)

    # x-axis: show in Mb
    mb_ticks = np.arange(
        np.ceil(start / 1e6) * 1e6,
        end + 1e6,
        1e6 if span_mb > 1.5 else 0.5e6
    )
    ax.set_xticks(mb_ticks)
    ax.set_xticklabels([f'{t/1e6:.1f}' for t in mb_ticks], fontsize=9)
    ax.set_xlabel(f'Genomic position — {region["chrom"]} (Mb)', fontsize=11)

    # y-axis: show in Mb (distance / 2)
    y_ticks = np.linspace(0, y_max, 5)
    ax.set_yticks(y_ticks)
    ax.set_yticklabels([f'{t/1e6:.2f}' for t in y_ticks], fontsize=8)
    ax.set_ylabel('Genomic distance / 2 (Mb)', fontsize=10)

    ax.set_title(
        f'{region["name"]}\n'
        f'{region["chrom"]}:{start:,}–{end:,}  '
        f'({span_mb:.2f} Mb, {n} bins @ {RES//1000} kb)  '
        f'|  Insulator: {region["ins_start"]:,}–{region["ins_end"]:,}',
        fontsize=11, fontweight='bold', pad=10
    )
    ax.legend(fontsize=9, loc='upper right')

    plt.colorbar(im, ax=ax, label='Contact frequency (log scale)',
                 shrink=0.7, pad=0.01)

    plt.tight_layout()
    path = f'media/{region["safe"]}_triangle.png'
    plt.savefig(path, dpi=220, bbox_inches='tight')
    print(f'  Saved: {path}')
    plt.close()
    return path


# ── Combined two-panel figure ─────────────────────────────────────────────────
def plot_combined(regions_data):
    """Both TAD regions as triangle plots stacked vertically for easy comparison."""
    fig, axes = plt.subplots(2, 1, figsize=(16, 11),
                             gridspec_kw={'hspace': 0.55})

    for ax, (region, matrix) in zip(axes, regions_data):
        n   = matrix.shape[0]
        start, end = region['start'], region['end']
        span_mb    = (end - start) / 1e6
        y_max      = (end - start) / 2

        tri    = rotate45(matrix)
        masked = np.ma.masked_invalid(tri)

        vals = matrix[~np.isnan(matrix) & (matrix > 0)]
        vmin = float(np.percentile(vals, 5))  if len(vals) else 1e-4
        vmax = float(np.percentile(vals, 99)) if len(vals) else 1.0

        cmap = plt.get_cmap('YlOrRd').copy()
        cmap.set_bad('white')
        norm = mcolors.LogNorm(vmin=max(vmin, 1e-5), vmax=vmax)

        extent = [start, end, 0, y_max]
        im = ax.imshow(masked, cmap=cmap, norm=norm, aspect='auto',
                       extent=extent, origin='lower', interpolation='nearest')

        ins_mid = (region['ins_start'] + region['ins_end']) / 2
        ax.axvline(ins_mid, color='#00bfff', lw=2.5, ls='--',
                   label='Insulator')
        ax.text(ins_mid, y_max * 0.96, region['ins_label'],
                ha='center', va='top', fontsize=8, color='#007acc',
                fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='#007acc',
                          alpha=0.85, lw=1))

        for pos in range(int(start / 500_000) * 500_000, end + 500_000, 500_000):
            if start < pos < end:
                ax.axvline(pos, color='#ccc', lw=0.4, ls=':')

        ax.set_xlim(start, end)
        ax.set_ylim(0, y_max)

        mb_ticks = np.arange(
            np.ceil(start / 1e6) * 1e6,
            end + 1e6,
            1e6 if span_mb > 2 else 0.5e6
        )
        ax.set_xticks(mb_ticks)
        ax.set_xticklabels([f'{t/1e6:.1f}' for t in mb_ticks], fontsize=9)
        ax.set_xlabel(f'{region["chrom"]} (Mb)', fontsize=10)

        y_ticks = np.linspace(0, y_max, 5)
        ax.set_yticks(y_ticks)
        ax.set_yticklabels([f'{t/1e6:.2f}' for t in y_ticks], fontsize=8)
        ax.set_ylabel('Distance/2 (Mb)', fontsize=9)

        ax.set_title(
            f'{region["name"]}   |   '
            f'{region["chrom"]}:{start:,}–{end:,}  ({span_mb:.2f} Mb)',
            fontsize=10, fontweight='bold'
        )
        ax.legend(fontsize=8, loc='upper right')
        plt.colorbar(im, ax=ax, label='Contact freq.', shrink=0.6, pad=0.01)

    fig.suptitle(
        'Full 2-TAD Pair Triangle Contact Maps — Mouse Micro-C (mm10, 5 kb)',
        fontsize=13, fontweight='bold', y=1.01
    )

    path = 'media/full_tad_triangles_combined.png'
    plt.savefig(path, dpi=220, bbox_inches='tight')
    print(f'  Saved: {path}')
    plt.close()
    return path


# ── Main ─────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    import os; os.makedirs('media', exist_ok=True)

    regions_data = []
    for region in REGIONS:
        print(f'\n{region["name"]}')
        matrix = clr.matrix(balance=True).fetch(region['query'])
        print(f'  Matrix: {matrix.shape[0]}×{matrix.shape[1]} '
              f'({np.sum(~np.isnan(matrix)):,} non-NaN values)')
        plot_triangle(region, matrix)
        regions_data.append((region, matrix))

    plot_combined(regions_data)
    print('\nDone.')
