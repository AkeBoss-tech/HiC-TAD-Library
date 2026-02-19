#!/usr/bin/env python3
"""
Deletion sensitivity scan — configurable region and deletion sizes.

Usage:
    python deletion_scan.py edward               # original insulator size
    python deletion_scan.py jingyun              # original insulator size
    python deletion_scan.py edward  --del-sizes 10 40 80   # 10/40/80 kb
    python deletion_scan.py jingyun --del-sizes 10 40 80

For each deletion size, 12 evenly-spaced deletions are simulated across the
1 Mb window.  The wild-type is predicted once and reused across all sizes.
Three metrics compare each deletion to WT:
  global_abs    — mean |Δ contact| (total reorganisation)
  cross_gain    — cross-boundary contact change (positive = TADs merging)
  ins_weakening — insulation score change at deletion site

Confirmed insulator coordinates (mm10):
  Edward:  chr12:27,333,532–27,336,455  (~3 kb)
  Jingyun: chr13:83,739,797–83,745,138  (~5.3 kb)

Cell type: EFO:0004038 (Mouse ESC)
"""

import os
import sys
import argparse
import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as mgs
from dotenv import load_dotenv

from alphagenome.data import genome
from alphagenome.models import dna_client

# ── Region configurations ──────────────────────────────────────────────────────
REGIONS = {
    'edward': dict(
        chrom        = 'chr12',
        del_start    = 27_333_532,
        del_end      = 27_336_455,
        name         = 'Edward',
        region_label = 'Edward chr12',
        safe         = 'edward_chr12',
    ),
    'jingyun': dict(
        chrom        = 'chr13',
        del_start    = 83_739_797,
        del_end      = 83_745_138,
        name         = 'Jingyun',
        region_label = 'Jingyun chr13',
        safe         = 'jingyun_chr13',
    ),
}

# ── Argument parsing ───────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description='Deletion sensitivity scan')
parser.add_argument('region', nargs='?', default='edward',
                    choices=list(REGIONS.keys()),
                    help='Region to scan (default: edward)')
parser.add_argument('--del-sizes', nargs='+', type=int, default=None,
                    metavar='KB',
                    help='Deletion sizes in kb, e.g. --del-sizes 10 40 80. '
                         'Default: use the actual insulator size.')
args = parser.parse_args()

cfg          = REGIONS[args.region]
CHROM        = cfg['chrom']
DEL_START    = cfg['del_start']
DEL_END      = cfg['del_end']
INSULATOR_SIZE = DEL_END - DEL_START   # actual biological insulator size
SAFE         = cfg['safe']
NAME         = cfg['name']
REGION_LABEL = cfg['region_label']

# Deletion sizes to test (in bp)
if args.del_sizes:
    DEL_SIZES_BP = [s * 1000 for s in args.del_sizes]
else:
    DEL_SIZES_BP = [INSULATOR_SIZE]   # original behaviour

# ── Model setup ───────────────────────────────────────────────────────────────
load_dotenv()
dna_model = dna_client.create(os.getenv('ALPHA_GENOME_API_KEY'))
print(f'Model initialized. Region: {REGION_LABEL}')
print(f'Deletion sizes: {[s//1000 for s in DEL_SIZES_BP]} kb')

ORGANISM  = dna_client.Organism.MUS_MUSCULUS
CELL_TYPE = 'EFO:0004038'
CELL_NAME = 'Mouse ESC'

# 1 Mb analysis window centred on the insulator
_mid      = (DEL_START + DEL_END) // 2
WIN_START = max(0, _mid - 2**19)
WIN_END   = _mid + 2**19
INTERVAL  = genome.Interval(chromosome=CHROM, start=WIN_START, end=WIN_END)

N_SITES = 12


# ── Helpers ───────────────────────────────────────────────────────────────────
def _rotate45(matrix):
    """Vectorised 45° rotation for classic triangle TAD plot."""
    n = matrix.shape[0]
    out = np.full((n, 2 * n), np.nan)
    rows, cols = np.triu_indices(n)
    out[cols - rows, rows + cols] = matrix[rows, cols]
    return out


def compute_metrics(wt, dm, del_bin, n, res):
    """
    Three complementary impact scores for one deletion.
    """
    diff = dm - wt
    global_abs = float(np.mean(np.abs(diff)))

    lo, hi = max(0, del_bin - 1), min(n, del_bin + 1)
    cross_gain = float(np.mean(dm[:lo, hi:]) - np.mean(wt[:lo, hi:]))

    w = max(2, int(50_000 / res))
    a, b = max(0, del_bin - w), del_bin
    c, d = del_bin, min(n, del_bin + w)
    ins_wt = float(np.mean(wt[a:b, c:d])) if b > a and d > c else 0.
    ins_dm = float(np.mean(dm[a:b, c:d])) if b > a and d > c else 0.
    ins_weakening = ins_dm - ins_wt

    return global_abs, cross_gain, ins_weakening


# ── Step 1: predict WT (once, shared across all deletion sizes) ───────────────
print(f'\nPredicting WT for {CHROM}:{WIN_START:,}-{WIN_END:,} ...')
wt_out = dna_model.predict_interval(
    interval=INTERVAL,
    requested_outputs={dna_client.OutputType.CONTACT_MAPS},
    ontology_terms=[CELL_TYPE],
    organism=ORGANISM,
)
wt  = wt_out.contact_maps.values[:, :, 0]
n   = wt.shape[0]
res = (WIN_END - WIN_START) / n
print(f'WT: {n}×{n} bins, resolution ≈ {res:,.0f} bp/bin\n')

actual_centre = (DEL_START + DEL_END) // 2


# ── Step 2–3: scan for each deletion size ─────────────────────────────────────
def run_scan(del_size_bp):
    """Run the 12-position scan for one deletion size. Returns list of record dicts."""
    print(f'\n── Scanning with {del_size_bp//1000} kb deletions ({del_size_bp:,} bp) ──')
    margin  = del_size_bp + int(res)
    centres = np.linspace(WIN_START + margin, WIN_END - margin - del_size_bp,
                          N_SITES, dtype=int)

    nearest = int(np.argmin(np.abs(centres - actual_centre)))
    centres[nearest] = actual_centre
    is_actual = np.zeros(N_SITES, dtype=bool)
    is_actual[nearest] = True

    records = []
    for i, centre in enumerate(centres):
        del_s   = int(centre)
        del_e   = del_s + del_size_bp
        del_bin = int((centre - WIN_START) / res)
        tag     = f'  ← {NAME.upper()} (actual insulator)' if is_actual[i] else ''
        print(f'  [{i+1:2d}/{N_SITES}]  {CHROM}:{del_s:,}–{del_e:,}{tag}')

        variant = genome.Variant(
            chromosome=CHROM,
            position=del_s - 1,
            reference_bases='N' * del_size_bp,
            alternate_bases='N',
        )
        del_out = dna_model.predict_variant(
            interval=INTERVAL,
            variant=variant,
            requested_outputs={dna_client.OutputType.CONTACT_MAPS},
            ontology_terms=[CELL_TYPE],
            organism=ORGANISM,
        )
        dm = del_out.alternate.contact_maps.values[:, :, 0]
        g, c, ins = compute_metrics(wt, dm, del_bin, n, res)
        records.append(dict(
            idx=i, centre=centre, del_s=del_s, del_e=del_e,
            del_bin=del_bin, dm=dm,
            global_abs=g, cross_gain=c, ins_weakening=ins,
            is_actual=bool(is_actual[i]),
        ))
        print(f'          global_abs={g:.5f}  cross_gain={c:+.5f}  '
              f'ins_weakening={ins:+.5f}')

    records.sort(key=lambda r: r['global_abs'], reverse=True)
    actual = next(r for r in records if r['is_actual'])
    rank   = records.index(actual) + 1
    print(f"  → {NAME}'s insulator ranks #{rank}/{N_SITES} by global impact.")
    return records


# Run all sizes
all_scans = {}   # del_size_bp → records
for dsz in DEL_SIZES_BP:
    all_scans[dsz] = run_scan(dsz)

print('\nAll scans complete.\n')


# ── Plotting helpers ──────────────────────────────────────────────────────────
def _size_label(bp):
    """Human-readable deletion size label."""
    return f'{bp//1000} kb' if bp >= 1000 else f'{bp} bp'


# ── Figure A: sensitivity profile (per size) ─────────────────────────────────
def plot_sensitivity_profile(records, del_size_bp):
    colours  = ['#e74c3c' if r['is_actual'] else '#3498db' for r in records]
    x_labels = [f"{r['centre'] / 1_000_000:.3f} Mb" for r in records]
    x        = np.arange(len(records))
    sl       = _size_label(del_size_bp)

    metrics = [
        ('global_abs',    'Mean |Δ contact|',
         'Total contact reorganisation (higher = bigger effect)'),
        ('cross_gain',    'Cross-TAD contact gain',
         'Mean contact gain across deletion site (positive = domains merging)'),
        ('ins_weakening', 'Insulation weakening',
         'Insulation score change at deletion site (positive = boundary weakened)'),
    ]

    fig, axes = plt.subplots(3, 1, figsize=(14, 11), sharex=True,
                             gridspec_kw={'hspace': 0.45})
    fig.suptitle(
        f'Deletion Sensitivity Scan — {REGION_LABEL} ({CELL_NAME})\n'
        f'Each bar = a {del_size_bp:,}-bp ({sl}) deletion centred at that position\n'
        f"Red bar = {NAME}'s actual insulator",
        fontsize=12, fontweight='bold'
    )

    for ax, (key, ylabel, desc) in zip(axes, metrics):
        vals = [r[key] for r in records]
        bars = ax.bar(x, vals, color=colours, edgecolor='white', linewidth=0.5)
        for j, r in enumerate(records):
            if r['is_actual']:
                bars[j].set_edgecolor('black')
                bars[j].set_linewidth(2)
                ypos = vals[j]
                ytxt = max(vals) * 0.85 if max(vals) != 0 else 0.001
                ax.annotate(f"{NAME}'s\ninsulator",
                            xy=(j, ypos), xytext=(j + 0.5, ytxt),
                            fontsize=7.5, color='#c0392b', fontweight='bold',
                            arrowprops=dict(arrowstyle='->', color='#c0392b', lw=1.2))
        ax.axhline(0, color='black', lw=0.6)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.set_title(desc, fontsize=8.5, color='#555', loc='left')
        ax.tick_params(axis='x', labelsize=7.5)

    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels(x_labels, rotation=35, ha='right')
    axes[-1].set_xlabel('Deletion centre (genomic position)', fontsize=10)

    path = f'media/deletion_scan_{SAFE}_{del_size_bp//1000}kb_sensitivity.png'
    plt.savefig(path, dpi=200, bbox_inches='tight')
    print(f'Saved: {path}')
    plt.close()
    return path


# ── Figure B: ranked summary per size ────────────────────────────────────────
def plot_ranked_summary(records, del_size_bp):
    sorted_by_pos  = sorted(records, key=lambda r: r['centre'])
    sorted_by_rank = sorted(records, key=lambda r: r['global_abs'], reverse=True)
    sl = _size_label(del_size_bp)

    cmap_pos = plt.get_cmap('coolwarm')
    pos_norm = plt.Normalize(WIN_START, WIN_END)

    fig = plt.figure(figsize=(16, 6))
    gs  = mgs.GridSpec(1, 2, width_ratios=[1.4, 1], wspace=0.35)

    ax_tri = fig.add_subplot(gs[0])
    start, end = WIN_START, WIN_END
    y_max = (end - start) / 2
    vmax  = np.percentile(wt[wt > 0], 99) if np.any(wt > 0) else 1
    cmap_obj = plt.get_cmap('Reds').copy(); cmap_obj.set_bad('white')
    ax_tri.imshow(np.ma.masked_invalid(_rotate45(wt)), cmap=cmap_obj,
                  aspect='auto', vmin=0, vmax=vmax, origin='lower',
                  extent=[start, end, 0, y_max], interpolation='nearest')
    for r in sorted_by_pos:
        col = '#e74c3c' if r['is_actual'] else cmap_pos(pos_norm(r['centre']))
        ax_tri.axvline(r['centre'] + del_size_bp // 2, color=col,
                       lw=2 if r['is_actual'] else 0.8,
                       ls='-' if r['is_actual'] else ':',
                       alpha=1.0 if r['is_actual'] else 0.6)
    ax_tri.set_xlim(start, end); ax_tri.set_ylim(0, y_max)
    ax_tri.set_title(f'Wild-type TAD structure\nScan positions for {sl} deletions '
                     f"(red = {NAME}'s insulator)", fontsize=9, fontweight='bold')
    ax_tri.set_xlabel('Genomic position (bp)', fontsize=8)
    ax_tri.set_ylabel('Distance/2 (bp)', fontsize=8)
    ax_tri.ticklabel_format(style='plain', axis='both')

    ax_bar = fig.add_subplot(gs[1])
    vals_ranked    = [r['global_abs'] for r in sorted_by_rank]
    colours_ranked = ['#e74c3c' if r['is_actual']
                      else cmap_pos(pos_norm(r['centre']))
                      for r in sorted_by_rank]
    y_pos = np.arange(len(records))
    bars  = ax_bar.barh(y_pos, vals_ranked, color=colours_ranked,
                        edgecolor='white', linewidth=0.5)
    for j, r in enumerate(sorted_by_rank):
        if r['is_actual']:
            bars[j].set_edgecolor('black'); bars[j].set_linewidth(2)
    ax_bar.set_yticks(y_pos)
    ax_bar.set_yticklabels(
        [f"{r['centre']/1_000_000:.3f} Mb" for r in sorted_by_rank], fontsize=7.5)
    ax_bar.invert_yaxis()
    ax_bar.set_xlabel('Mean |Δ contact| (impact score)', fontsize=9)
    ax_bar.set_title('Impact ranking\n(top = biggest structural change)',
                     fontsize=9, fontweight='bold')
    ax_bar.axvline(0, color='black', lw=0.5)

    fig.suptitle(
        f'Deletion Scan Summary — {CHROM} | {CELL_NAME} | {sl} deletions\n'
        f"Red = {NAME}'s actual insulator at {actual_centre / 1_000_000:.3f} Mb",
        fontsize=11, fontweight='bold'
    )
    path = f'media/deletion_scan_{SAFE}_{del_size_bp//1000}kb_ranked_summary.png'
    plt.savefig(path, dpi=200, bbox_inches='tight')
    print(f'Saved: {path}')
    plt.close()
    return path


# ── Figure C: triangle gallery per size ──────────────────────────────────────
def plot_triangle_gallery(records, del_size_bp):
    actual   = next(r for r in records if r['is_actual'])
    top2     = [r for r in records if not r['is_actual']][:2]
    others   = [r for r in records if not r['is_actual'] and r not in top2]
    step     = max(1, len(others) // 3)
    controls = others[::step][:3]
    selected = [actual] + top2 + controls

    start, end = WIN_START, WIN_END
    y_max  = (end - start) / 2
    extent = [start, end, 0, y_max]
    vmax   = np.percentile(wt[wt > 0], 99) if np.any(wt > 0) else 1
    sl     = _size_label(del_size_bp)

    ncols    = 3
    nrows    = int(np.ceil((len(selected) + 1) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 4.5 * nrows),
                             gridspec_kw={'hspace': 0.55, 'wspace': 0.35})
    axes_flat = axes.flatten()
    colour_map = {'actual': '#e74c3c', 'top': '#e67e22', 'ctrl': '#2980b9'}

    fig.suptitle(
        f'Triangle TAD Views — {sl} Deletion Scan ({CELL_NAME})\n'
        f"Red border = {NAME}'s actual insulator  |  "
        f'Orange = highest-impact non-insulator  |  Blue = interior controls',
        fontsize=11, fontweight='bold'
    )

    # WT panel
    cmap_wt = plt.get_cmap('Reds').copy(); cmap_wt.set_bad('white')
    im0 = axes_flat[0].imshow(np.ma.masked_invalid(_rotate45(wt)), cmap=cmap_wt,
                               aspect='auto', vmin=0, vmax=vmax,
                               extent=extent, origin='lower', interpolation='nearest')
    axes_flat[0].set_xlim(start, end); axes_flat[0].set_ylim(0, y_max)
    axes_flat[0].set_title(f'Wild-type (WT)\n{CHROM}:{WIN_START:,}–{WIN_END:,}',
                            fontsize=8.5, fontweight='bold', pad=3)
    axes_flat[0].set_xlabel('Genomic position (bp)', fontsize=7.5)
    axes_flat[0].set_ylabel('Distance/2 (bp)', fontsize=7.5)
    axes_flat[0].ticklabel_format(style='plain', axis='both')
    plt.colorbar(im0, ax=axes_flat[0], shrink=0.5, label='Contact freq.')
    for spine in axes_flat[0].spines.values():
        spine.set_edgecolor('black'); spine.set_linewidth(3)

    for idx, r in enumerate(selected, start=1):
        if idx >= len(axes_flat): break
        ax   = axes_flat[idx]
        diff = r['dm'] - wt
        dvlim = np.percentile(np.abs(diff), 99)
        if r['is_actual']:
            border, role = colour_map['actual'], "ACTUAL INSULATOR"
        elif r in top2:
            border, role = colour_map['top'], f"High impact #{top2.index(r)+1}"
        else:
            border, role = colour_map['ctrl'], "Interior control"
        label = (f"{role}\nDel: {r['del_s']:,}–{r['del_e']:,}\n"
                 f"Δ={r['global_abs']:.4f}  cross={r['cross_gain']:+.4f}")
        cmap_d = plt.get_cmap('RdBu_r').copy(); cmap_d.set_bad('white')
        im = ax.imshow(np.ma.masked_invalid(_rotate45(diff)), cmap=cmap_d,
                       aspect='auto', vmin=-dvlim, vmax=dvlim,
                       extent=extent, origin='lower', interpolation='nearest')
        ax.axvline(r['del_s'], color='cyan', lw=1.2, ls='--', alpha=0.9)
        ax.axvline(r['del_e'], color='cyan', lw=1.2, ls='--', alpha=0.9)
        ax.set_xlim(start, end); ax.set_ylim(0, y_max)
        ax.set_title(label, fontsize=7.5, fontweight='bold', pad=3)
        ax.set_xlabel('Genomic position (bp)', fontsize=7.5)
        ax.set_ylabel('Distance/2 (bp)', fontsize=7.5)
        ax.ticklabel_format(style='plain', axis='both')
        plt.colorbar(im, ax=ax, shrink=0.5, label='Δ Contact freq.')
        for spine in ax.spines.values():
            spine.set_edgecolor(border); spine.set_linewidth(3)

    for ax in axes_flat[len(selected) + 1:]:
        ax.set_visible(False)

    path = f'media/deletion_scan_{SAFE}_{del_size_bp//1000}kb_triangle_gallery.png'
    plt.savefig(path, dpi=200, bbox_inches='tight')
    print(f'Saved: {path}')
    plt.close()
    return path


# ── Figure D: difference-map montage per size ─────────────────────────────────
def plot_diff_montage(records, del_size_bp):
    ncols = 4
    nrows = int(np.ceil(N_SITES / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.5 * nrows),
                             gridspec_kw={'hspace': 0.6, 'wspace': 0.35})
    axes_flat = axes.flatten()

    start, end = WIN_START, WIN_END
    y_max  = (end - start) / 2
    extent = [start, end, 0, y_max]
    all_diffs   = np.concatenate([np.abs(r['dm'] - wt).ravel() for r in records])
    global_vlim = float(np.percentile(all_diffs, 99))
    sorted_records = sorted(records, key=lambda r: r['centre'])
    ranked = sorted(records, key=lambda r: r['global_abs'], reverse=True)

    for k, r in enumerate(sorted_records):
        ax   = axes_flat[k]
        diff = r['dm'] - wt
        cmap_obj = plt.get_cmap('RdBu_r').copy(); cmap_obj.set_bad('white')
        ax.imshow(np.ma.masked_invalid(_rotate45(diff)), cmap=cmap_obj, aspect='auto',
                  vmin=-global_vlim, vmax=global_vlim,
                  extent=extent, origin='lower', interpolation='nearest')
        ax.axvline(r['del_s'], color='cyan', lw=1, ls='--', alpha=0.9)
        ax.axvline(r['del_e'], color='cyan', lw=1, ls='--', alpha=0.9)
        ax.set_xlim(start, end); ax.set_ylim(0, y_max)
        rank_n = ranked.index(r) + 1
        border = '#e74c3c' if r['is_actual'] else '#888'
        title = (f"Rank #{rank_n}{'  ← ACTUAL' if r['is_actual'] else ''}\n"
                 f"{r['centre'] / 1_000_000:.3f} Mb  Δ={r['global_abs']:.4f}")
        ax.set_title(title, fontsize=7,
                     fontweight='bold' if r['is_actual'] else 'normal')
        ax.tick_params(labelsize=5)
        ax.ticklabel_format(style='plain', axis='both')
        for spine in ax.spines.values():
            spine.set_edgecolor(border)
            spine.set_linewidth(2.5 if r['is_actual'] else 0.5)

    for ax in axes_flat[len(sorted_records):]:
        ax.set_visible(False)

    sl = _size_label(del_size_bp)
    actual_label = f"{NAME}'s actual insulator"
    fig.suptitle(
        f'All {N_SITES} Deletion Sites — Difference Maps (del − WT) | {sl}\n'
        f'{CELL_NAME}  |  Colour scale ±{global_vlim:.3f}  |  '
        f'Red border = {actual_label}',
        fontsize=11, fontweight='bold'
    )
    path = f'media/deletion_scan_{SAFE}_{del_size_bp//1000}kb_montage.png'
    plt.savefig(path, dpi=200, bbox_inches='tight')
    print(f'Saved: {path}')
    plt.close()
    return path


# ── Figure E: cross-size comparison (only when multiple sizes) ────────────────
def plot_cross_size_comparison(all_scans):
    """
    Grid showing how global_abs changes with deletion size.
    Rows = deletion sizes; columns = genomic positions (sorted).
    The actual insulator column is highlighted in red.
    Also shows the rank of the actual insulator at each size.
    """
    sizes     = sorted(all_scans.keys())
    if len(sizes) < 2:
        return None

    # Get positions from first scan (same for all sizes by design)
    first = sorted(all_scans[sizes[0]], key=lambda r: r['centre'])
    positions_mb = [r['centre'] / 1_000_000 for r in first]
    actual_idx = next(i for i, r in enumerate(first) if r['is_actual'])

    # Build matrix: rows=sizes, cols=positions (sorted genomically)
    mat_abs = np.zeros((len(sizes), len(first)))
    mat_ins = np.zeros((len(sizes), len(first)))
    ranks   = []   # rank of actual insulator at each deletion size

    for si, sz in enumerate(sizes):
        records_sorted = sorted(all_scans[sz], key=lambda r: r['centre'])
        records_ranked = sorted(all_scans[sz], key=lambda r: r['global_abs'], reverse=True)
        for pi, r in enumerate(records_sorted):
            mat_abs[si, pi] = r['global_abs']
            mat_ins[si, pi] = r['ins_weakening']
        actual_r = next(r for r in records_ranked if r['is_actual'])
        ranks.append(records_ranked.index(actual_r) + 1)

    fig, axes = plt.subplots(3, 1, figsize=(14, 12),
                             gridspec_kw={'hspace': 0.55})

    # Panel 1: heatmap of global_abs
    ax = axes[0]
    im = ax.imshow(mat_abs, aspect='auto', cmap='YlOrRd', origin='upper')
    ax.set_xticks(range(len(positions_mb)))
    ax.set_xticklabels([f'{p:.3f}' for p in positions_mb], rotation=40, ha='right', fontsize=7)
    ax.set_yticks(range(len(sizes)))
    ax.set_yticklabels([_size_label(s) for s in sizes], fontsize=9)
    ax.set_xlabel('Deletion centre (Mb)', fontsize=9)
    ax.set_ylabel('Deletion size', fontsize=9)
    ax.set_title('Mean |Δ contact| at each position × size\n'
                 '(higher = more structural disruption)', fontsize=9, loc='left')
    plt.colorbar(im, ax=ax, label='Mean |Δ contact|', shrink=0.6)
    # Highlight actual insulator column
    ax.axvline(actual_idx - 0.5, color='#e74c3c', lw=2)
    ax.axvline(actual_idx + 0.5, color='#e74c3c', lw=2)
    ax.text(actual_idx, -0.7, f"{NAME}'s\ninsulator", ha='center', va='top',
            fontsize=7, color='#e74c3c', fontweight='bold', transform=ax.get_xaxis_transform())

    # Panel 2: heatmap of insulation weakening
    ax2 = axes[1]
    vlim = np.percentile(np.abs(mat_ins), 98)
    im2  = ax2.imshow(mat_ins, aspect='auto', cmap='RdBu_r',
                      vmin=-vlim, vmax=vlim, origin='upper')
    ax2.set_xticks(range(len(positions_mb)))
    ax2.set_xticklabels([f'{p:.3f}' for p in positions_mb], rotation=40, ha='right', fontsize=7)
    ax2.set_yticks(range(len(sizes)))
    ax2.set_yticklabels([_size_label(s) for s in sizes], fontsize=9)
    ax2.set_xlabel('Deletion centre (Mb)', fontsize=9)
    ax2.set_ylabel('Deletion size', fontsize=9)
    ax2.set_title('Insulation weakening at each position × size\n'
                  '(positive = boundary lost; negative = boundary strengthened)',
                  fontsize=9, loc='left')
    plt.colorbar(im2, ax=ax2, label='Insulation weakening', shrink=0.6)
    ax2.axvline(actual_idx - 0.5, color='#e74c3c', lw=2)
    ax2.axvline(actual_idx + 0.5, color='#e74c3c', lw=2)
    ax2.text(actual_idx, -0.7, f"{NAME}'s\ninsulator", ha='center', va='top',
             fontsize=7, color='#e74c3c', fontweight='bold',
             transform=ax2.get_xaxis_transform())

    # Panel 3: rank of actual insulator across sizes
    ax3 = axes[2]
    ax3.plot([_size_label(s) for s in sizes], ranks, 'o-',
             color='#e74c3c', lw=2.5, ms=8, markerfacecolor='white',
             markeredgewidth=2)
    for xi, (s, rk) in enumerate(zip(sizes, ranks)):
        ax3.annotate(f'#{rk}', (xi, rk), textcoords='offset points',
                     xytext=(6, 4), fontsize=9, color='#c0392b', fontweight='bold')
    ax3.set_ylim(0, N_SITES + 1)
    ax3.invert_yaxis()
    ax3.axhline(1, color='#aaa', lw=0.8, ls='--')
    ax3.set_xlabel('Deletion size', fontsize=10)
    ax3.set_ylabel('Rank of actual insulator\n(1 = biggest global impact)', fontsize=9)
    ax3.set_title(f"{NAME}'s insulator rank vs. deletion size\n"
                  f"(lower = insulator causes more global disruption than controls)",
                  fontsize=9, loc='left')
    ax3.grid(axis='y', lw=0.4, alpha=0.5)

    fig.suptitle(
        f'Cross-Size Comparison — {REGION_LABEL} | {CELL_NAME}\n'
        f'How deletion size affects which genomic positions cause the most disruption',
        fontsize=12, fontweight='bold'
    )
    path = f'media/deletion_scan_{SAFE}_cross_size_comparison.png'
    plt.savefig(path, dpi=200, bbox_inches='tight')
    print(f'Saved: {path}')
    plt.close()
    return path


# ── HTML report ────────────────────────────────────────────────────────────────
def generate_scan_html(size_paths, all_scans, cross_path=None):
    """Generate combined HTML report for all deletion sizes."""
    now    = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
    sizes  = sorted(all_scans.keys())
    name_poss = NAME + "'s"

    def img(src, cap='', w='100%'):
        return (f'<figure style="margin:0 0 16px 0">'
                f'<img src="{src}" style="width:{w};border:1px solid #ddd;'
                f'border-radius:6px">'
                f'<figcaption style="font-size:.82em;color:#555;margin-top:4px">'
                f'{cap}</figcaption></figure>')

    # Per-size tables
    size_sections = ''
    for sz in sizes:
        sl      = _size_label(sz)
        records = all_scans[sz]
        ranked  = sorted(records, key=lambda r: r['global_abs'], reverse=True)
        rows    = ''
        for rank, r in enumerate(ranked, 1):
            flag  = ' ⭐ ACTUAL' if r['is_actual'] else ''
            style = 'background:#fdecea;font-weight:bold' if r['is_actual'] else ''
            rows += (f'<tr style="{style}">'
                     f'<td>#{rank}</td><td>{r["centre"]:,}</td>'
                     f'<td>{r["del_s"]:,}–{r["del_e"]:,}</td>'
                     f'<td>{r["global_abs"]:.5f}{flag}</td>'
                     f'<td>{r["cross_gain"]:+.5f}</td>'
                     f'<td>{r["ins_weakening"]:+.5f}</td></tr>\n')

        p = size_paths[sz]
        cap_sum  = f'WT triangle with {sl} scan positions. Right: impact ranking.'
        cap_sens = (f'Three metrics at each position for {sl} deletions. '
                    f'Red bar = {name_poss} actual insulator.')
        cap_gal  = f'Triangle difference maps for selected sites ({sl} deletions).'
        cap_mon  = f'All {N_SITES} sites — {sl} deletions — shared colour scale.'

        size_sections += f"""
<h2>{sl} Deletions</h2>
<table>
  <tr><th>Rank</th><th>Centre (bp)</th><th>Deletion range</th>
      <th>Mean |&Delta; contact|</th><th>Cross-TAD gain</th><th>Ins. weakening</th></tr>
  {rows}
</table>
{img(p['summary'],     cap_sum)}
{img(p['sensitivity'], cap_sens)}
{img(p['gallery'],     cap_gal)}
{img(p['montage'],     cap_mon)}
"""

    cross_html = ''
    if cross_path:
        cross_html = f"""
<h2>Cross-Size Comparison</h2>
<p style="font-size:.88em;color:#444">
  How deletion size affects the spatial sensitivity pattern.
  Top: global impact heatmap (position × size).
  Middle: insulation weakening heatmap.
  Bottom: rank of {name_poss} actual insulator at each deletion size.
</p>
{img(cross_path, f'Cross-size comparison for {REGION_LABEL}.')}
"""

    size_list  = ', '.join(_size_label(s) for s in sizes)
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>Deletion Scan — {NAME} {CHROM}</title>
  <style>
    body{{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Arial,sans-serif;
         max-width:1300px;margin:0 auto;padding:24px;background:#f0f2f5;color:#2c3e50}}
    h1{{font-size:1.6em;border-bottom:3px solid #2980b9;padding-bottom:8px}}
    h2{{font-size:1.2em;color:#1a5276;margin:28px 0 8px}}
    .meta{{font-size:.85em;color:#888;margin-bottom:20px}}
    .info{{background:#eaf4fb;padding:12px 16px;border-radius:6px;
           font-size:.88em;margin-bottom:18px;line-height:1.75}}
    table{{border-collapse:collapse;width:100%;font-size:.85em;margin-bottom:24px}}
    th{{background:#2980b9;color:white;padding:7px 10px;text-align:left}}
    td{{padding:6px 10px;border-bottom:1px solid #ddd}}
    tr:hover td{{background:#f5faff}}
    .legend{{background:#fdfefe;border:1px solid #d5d8dc;border-radius:8px;
             padding:14px;font-size:.85em;margin:16px 0}}
  </style>
</head>
<body>
<h1>Deletion Sensitivity Scan — {NAME} {CHROM}</h1>
<p class="meta">
  {CHROM} | {CELL_NAME} ({CELL_TYPE}) | mm10 |
  {N_SITES} deletion sites | Sizes tested: {size_list} | Generated: {now}
</p>

<div class="info">
  <b>Hypothesis:</b> Deletions <em>at</em> a TAD boundary cause the largest structural
  change (TADs merge); deletions <em>within</em> a TAD body cause little change.<br>
  <b>Method:</b> Predict WT once; for each deletion size and each of {N_SITES}
  evenly-spaced positions, simulate a deletion and compare to WT.<br>
  <b>{NAME}'s actual insulator:</b> {CHROM}:{DEL_START:,}–{DEL_END:,}
  ({INSULATOR_SIZE:,} bp) — highlighted in red.
</div>

<div class="legend">
  <b>Metrics explained</b><br>
  &bull; <b>Mean |&Delta; contact|</b> — total reorganisation; higher = more change.<br>
  &bull; <b>Cross-TAD contact gain</b> — positive = domains merging.<br>
  &bull; <b>Insulation weakening</b> — positive = boundary lost.
</div>

{cross_html}
{size_sections}

</body>
</html>
"""
    suffix = '_'.join(f'{s//1000}kb' for s in sizes)
    out    = f'deletion_scan_{SAFE}_{suffix}_report.html'
    with open(out, 'w') as fh:
        fh.write(html)
    print(f'HTML report saved: {out}')
    return out


# ── Run all plots ─────────────────────────────────────────────────────────────
if __name__ == '__main__':
    os.makedirs('media', exist_ok=True)

    size_paths = {}
    for dsz, records in all_scans.items():
        size_paths[dsz] = {
            'sensitivity': plot_sensitivity_profile(records, dsz),
            'gallery':     plot_triangle_gallery(records, dsz),
            'montage':     plot_diff_montage(records, dsz),
            'summary':     plot_ranked_summary(records, dsz),
        }

    cross_path = plot_cross_size_comparison(all_scans)

    html_path = generate_scan_html(size_paths, all_scans, cross_path)

    # Print summary table
    print('\n── Results by size ──────────────────────────────────────────────────')
    for dsz in sorted(all_scans.keys()):
        records  = all_scans[dsz]
        ranked   = sorted(records, key=lambda r: r['global_abs'], reverse=True)
        actual   = next(r for r in ranked if r['is_actual'])
        rank     = ranked.index(actual) + 1
        ins_wk   = actual['ins_weakening']
        print(f'  {_size_label(dsz):>6}  →  rank #{rank}/{N_SITES}  '
              f'global_abs={actual["global_abs"]:.5f}  '
              f'ins_weakening={ins_wk:+.5f}')

    print(f'\nHTML report: {html_path}')
    print('Figures saved to media/.')
