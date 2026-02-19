#!/usr/bin/env python3
"""
Deletion sensitivity scan across Edward's chr12 region.

Tests N evenly-spaced deletions — all the same size as Edward's real insulator
(~3 kb) — across the 1 Mb analysis window.  Compares each one against the
wild-type contact map to reveal which genomic positions, when deleted, produce
the largest structural change.

The hypothesis:
  Deletions AT a TAD boundary  → large contact reorganisation (TADs merge)
  Deletions WITHIN a TAD body → small change (domain interior is redundant)

Edward's confirmed coordinates (mm10):
  Insulator: chr12:27,333,532–27,336,455  (~3 kb)
  Full 2-TAD span: chr12:26,440,002–28,560,000  (2.1 Mb)

Cell type: EFO:0004038 (Mouse ESC) — showed the most loop structure for chr12.
"""

import os
import datetime
import textwrap
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as mgs
from matplotlib.patches import FancyArrowPatch
from dotenv import load_dotenv

from alphagenome.data import genome
from alphagenome.models import dna_client

# ── Setup ─────────────────────────────────────────────────────────────────────
load_dotenv()
dna_model = dna_client.create(os.getenv('ALPHA_GENOME_API_KEY'))
print('Model initialized.')

ORGANISM  = dna_client.Organism.MUS_MUSCULUS
CELL_TYPE = 'EFO:0004038'
CELL_NAME = 'Mouse ESC'

# Edward's insulator (mm10)
ED_CHROM     = 'chr12'
ED_DEL_START = 27_333_532
ED_DEL_END   = 27_336_455
DEL_SIZE     = ED_DEL_END - ED_DEL_START   # 2,923 bp

# 1 Mb analysis window centred on Edward's deletion
_mid      = (ED_DEL_START + ED_DEL_END) // 2
WIN_START = max(0, _mid - 2**19)
WIN_END   = _mid + 2**19
INTERVAL  = genome.Interval(chromosome=ED_CHROM, start=WIN_START, end=WIN_END)

# Number of candidate deletion centres to test
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

    global_abs    : mean |Δ contact| across the entire map — total reorganisation
    cross_gain    : mean contact increase between genomic left/right of deletion —
                    positive = TADs merging
    ins_weakening : insulation score change at the deletion site —
                    positive = boundary weakened (less insulation)
    """
    diff = dm - wt

    global_abs = float(np.mean(np.abs(diff)))

    # Cross-boundary: all contacts between genomic left and right of deletion site
    lo, hi = max(0, del_bin - 1), min(n, del_bin + 1)
    cross_gain = float(np.mean(dm[:lo, hi:]) - np.mean(wt[:lo, hi:]))

    # Insulation score: mean contact in a diamond centred on del_bin
    w = max(2, int(50_000 / res))   # 50 kb window
    a, b = max(0, del_bin - w), del_bin
    c, d = del_bin, min(n, del_bin + w)
    ins_wt = float(np.mean(wt[a:b, c:d])) if b > a and d > c else 0.
    ins_dm = float(np.mean(dm[a:b, c:d])) if b > a and d > c else 0.
    ins_weakening = ins_dm - ins_wt

    return global_abs, cross_gain, ins_weakening


# ── Step 1: predict WT ────────────────────────────────────────────────────────
print(f'\nPredicting WT for {ED_CHROM}:{WIN_START:,}-{WIN_END:,} ...')
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


# ── Step 2: build scan positions ──────────────────────────────────────────────
margin  = DEL_SIZE + int(res)
centres = np.linspace(WIN_START + margin, WIN_END - margin - DEL_SIZE,
                      N_SITES, dtype=int)

# Replace the nearest site with Edward's actual centre
ed_centre  = (ED_DEL_START + ED_DEL_END) // 2
nearest    = int(np.argmin(np.abs(centres - ed_centre)))
centres[nearest] = ed_centre
is_actual  = np.zeros(N_SITES, dtype=bool)
is_actual[nearest] = True


# ── Step 3: scan ──────────────────────────────────────────────────────────────
records = []
for i, centre in enumerate(centres):
    del_s   = int(centre)
    del_e   = del_s + DEL_SIZE
    del_bin = int((centre - WIN_START) / res)
    tag     = '  ← EDWARD (actual insulator)' if is_actual[i] else ''
    print(f'  [{i+1:2d}/{N_SITES}]  {ED_CHROM}:{del_s:,}–{del_e:,}{tag}')

    variant = genome.Variant(
        chromosome=ED_CHROM,
        position=del_s - 1,
        reference_bases='N' * DEL_SIZE,
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

print('\nScan complete.\n')

# Rank by global_abs
records.sort(key=lambda r: r['global_abs'], reverse=True)
rank_labels = [f"#{k+1}" for k in range(len(records))]


# ── Figure 1: sensitivity profile ────────────────────────────────────────────
def plot_sensitivity_profile(records):
    """Bar charts of all three metrics at each deletion site."""
    positions = [r['centre'] for r in records]
    colours   = ['#e74c3c' if r['is_actual'] else '#3498db' for r in records]
    x_labels  = [f"{r['centre'] // 1_000_000:.3f} Mb" for r in records]
    x         = np.arange(len(records))

    metrics = [
        ('global_abs',    'Mean |Δ contact|',       'Total contact reorganisation\n(higher = bigger effect)'),
        ('cross_gain',    'Cross-TAD contact gain',  'Mean contact gain across deletion site\n(positive = domains merging)'),
        ('ins_weakening', 'Insulation weakening',    'Insulation score change at deletion site\n(positive = boundary weakened)'),
    ]

    fig, axes = plt.subplots(3, 1, figsize=(14, 11), sharex=True,
                             gridspec_kw={'hspace': 0.45})
    fig.suptitle(
        f'Deletion Sensitivity Scan — Edward chr12 ({CELL_NAME})\n'
        f'Each bar = a {DEL_SIZE:,}-bp deletion centred at that genomic position\n'
        f'Red bar = Edward\'s actual insulator',
        fontsize=12, fontweight='bold'
    )

    for ax, (key, ylabel, desc) in zip(axes, metrics):
        vals = [r[key] for r in records]
        bars = ax.bar(x, vals, color=colours, edgecolor='white', linewidth=0.5)
        # Mark the actual insulator
        for j, r in enumerate(records):
            if r['is_actual']:
                bars[j].set_edgecolor('black')
                bars[j].set_linewidth(2)
                ax.annotate("Edward's\ninsulator", xy=(j, vals[j]),
                            xytext=(j + 0.5, max(vals) * 0.85),
                            fontsize=7.5, color='#c0392b', fontweight='bold',
                            arrowprops=dict(arrowstyle='->', color='#c0392b', lw=1.2))
        ax.axhline(0, color='black', lw=0.6)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.set_title(desc, fontsize=8.5, color='#555', loc='left')
        ax.tick_params(axis='x', labelsize=7.5)

    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels(x_labels, rotation=35, ha='right')
    axes[-1].set_xlabel('Deletion centre (genomic position)', fontsize=10)

    path = 'media/deletion_scan_sensitivity.png'
    plt.savefig(path, dpi=200, bbox_inches='tight')
    print(f'Saved: {path}')
    plt.close()
    return path


# ── Figure 2: triangle TAD gallery ───────────────────────────────────────────
def plot_triangle_gallery(records, wt):
    """
    Triangle TAD views for 6 selected deletion sites + WT.
    Selects: actual insulator, top-2 non-actual, 3 random interior sites.
    """
    actual  = next(r for r in records if r['is_actual'])
    # top-2 by global impact that are not the actual site
    top2    = [r for r in records if not r['is_actual']][:2]
    # 3 "interior" control sites: evenly pick from remainder, avoid ends
    others  = [r for r in records if not r['is_actual'] and r not in top2]
    step    = max(1, len(others) // 3)
    controls = others[::step][:3]

    selected = [actual] + top2 + controls   # up to 6 sites
    N = len(selected)

    start, end = WIN_START, WIN_END
    y_max  = (end - start) / 2
    extent = [start, end, 0, y_max]
    vmax   = np.percentile(wt[wt > 0], 99) if np.any(wt > 0) else 1

    ncols = 3
    nrows = int(np.ceil((N + 1) / ncols))   # +1 for WT panel

    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 4.5 * nrows),
                             gridspec_kw={'hspace': 0.55, 'wspace': 0.35})
    axes_flat = axes.flatten()

    fig.suptitle(
        f'Triangle TAD Views — Deletion Scan ({CELL_NAME})\n'
        f'Red border = Edward\'s actual insulator  |  '
        f'Orange = highest-impact non-insulator  |  Blue = interior controls',
        fontsize=11, fontweight='bold'
    )

    def _tri_panel(ax, mat, title, cmap='Reds', vmin=0, vmax_use=None, border=None):
        cmap_obj = plt.get_cmap(cmap).copy()
        cmap_obj.set_bad('white')
        masked = np.ma.masked_invalid(_rotate45(mat))
        im = ax.imshow(masked, cmap=cmap_obj, aspect='auto',
                       vmin=vmin, vmax=vmax_use or vmax,
                       extent=extent, origin='lower', interpolation='nearest')
        ax.set_xlim(start, end)
        ax.set_ylim(0, y_max)
        ax.set_title(title, fontsize=8.5, fontweight='bold', pad=3)
        ax.set_xlabel('Genomic position (bp)', fontsize=7.5)
        ax.set_ylabel('Distance/2 (bp)', fontsize=7.5)
        ax.ticklabel_format(style='plain', axis='both')
        plt.colorbar(im, ax=ax, shrink=0.5, label='Contact freq.')
        if border:
            for spine in ax.spines.values():
                spine.set_edgecolor(border)
                spine.set_linewidth(3)

    # WT first
    _tri_panel(axes_flat[0], wt,
               f'Wild-type (WT)\n{ED_CHROM}:{WIN_START:,}–{WIN_END:,}',
               border='black')

    colour_map = {
        'actual': '#e74c3c',
        'top':    '#e67e22',
        'ctrl':   '#2980b9',
    }
    panel_idx = 1
    for r in selected:
        if panel_idx >= len(axes_flat):
            break
        ax = axes_flat[panel_idx]
        diff = r['dm'] - wt
        dvlim = np.percentile(np.abs(diff), 99)

        if r['is_actual']:
            border, role = colour_map['actual'], "ACTUAL INSULATOR"
        elif r in top2:
            border, role = colour_map['top'], f"High impact #{top2.index(r)+1}"
        else:
            border, role = colour_map['ctrl'], "Interior control"

        label = (
            f"{role}\n"
            f"Del: {r['del_s']:,}–{r['del_e']:,}\n"
            f"Δ={r['global_abs']:.4f}  cross={r['cross_gain']:+.4f}"
        )
        # Show difference map (del − WT) for clearer visual
        cmap_obj = plt.get_cmap('RdBu_r').copy()
        cmap_obj.set_bad('white')
        masked = np.ma.masked_invalid(_rotate45(diff))
        im = ax.imshow(masked, cmap=cmap_obj, aspect='auto',
                       vmin=-dvlim, vmax=dvlim,
                       extent=extent, origin='lower', interpolation='nearest')
        ax.axvline(r['del_s'], color='cyan', lw=1.2, ls='--', alpha=0.9)
        ax.axvline(r['del_e'], color='cyan', lw=1.2, ls='--', alpha=0.9)
        ax.set_xlim(start, end)
        ax.set_ylim(0, y_max)
        ax.set_title(label, fontsize=7.5, fontweight='bold', pad=3)
        ax.set_xlabel('Genomic position (bp)', fontsize=7.5)
        ax.set_ylabel('Distance/2 (bp)', fontsize=7.5)
        ax.ticklabel_format(style='plain', axis='both')
        plt.colorbar(im, ax=ax, shrink=0.5, label='Δ Contact freq.')
        for spine in ax.spines.values():
            spine.set_edgecolor(border)
            spine.set_linewidth(3)

        panel_idx += 1

    # Hide unused panels
    for ax in axes_flat[panel_idx:]:
        ax.set_visible(False)

    path = 'media/deletion_scan_triangle_gallery.png'
    plt.savefig(path, dpi=200, bbox_inches='tight')
    print(f'Saved: {path}')
    plt.close()
    return path


# ── Figure 3: compact difference-map montage (all sites) ─────────────────────
def plot_diff_montage(records, wt):
    """
    Small difference-map thumbnail for every scanned site.
    Laid out in a grid so you can scan all at once.
    """
    ncols = 4
    nrows = int(np.ceil(N_SITES / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(5 * ncols, 4.5 * nrows),
                             gridspec_kw={'hspace': 0.6, 'wspace': 0.35})
    axes_flat = axes.flatten()

    start, end = WIN_START, WIN_END
    y_max  = (end - start) / 2
    extent = [start, end, 0, y_max]

    # Global colour scale so all panels are comparable
    all_diffs = np.concatenate([np.abs(r['dm'] - wt).ravel() for r in records])
    global_vlim = float(np.percentile(all_diffs, 99))

    # Sort back to genomic order for the montage
    sorted_records = sorted(records, key=lambda r: r['centre'])

    for k, r in enumerate(sorted_records):
        ax = axes_flat[k]
        diff = r['dm'] - wt
        cmap_obj = plt.get_cmap('RdBu_r').copy()
        cmap_obj.set_bad('white')
        masked = np.ma.masked_invalid(_rotate45(diff))
        ax.imshow(masked, cmap=cmap_obj, aspect='auto',
                  vmin=-global_vlim, vmax=global_vlim,
                  extent=extent, origin='lower', interpolation='nearest')
        ax.axvline(r['del_s'], color='cyan', lw=1, ls='--', alpha=0.9)
        ax.axvline(r['del_e'], color='cyan', lw=1, ls='--', alpha=0.9)
        ax.set_xlim(start, end)
        ax.set_ylim(0, y_max)
        # Rank by impact (1 = biggest)
        rank = sorted(records, key=lambda x: x['global_abs'], reverse=True).index(r) + 1
        border = '#e74c3c' if r['is_actual'] else '#888'
        title = (
            f"Impact rank #{rank}{'  ← ACTUAL' if r['is_actual'] else ''}\n"
            f"{r['centre'] // 1_000_000:.3f} Mb  Δ={r['global_abs']:.4f}"
        )
        ax.set_title(title, fontsize=7, fontweight='bold' if r['is_actual'] else 'normal')
        ax.tick_params(labelsize=5)
        ax.ticklabel_format(style='plain', axis='both')
        for spine in ax.spines.values():
            spine.set_edgecolor(border)
            spine.set_linewidth(2.5 if r['is_actual'] else 0.5)

    for ax in axes_flat[len(sorted_records):]:
        ax.set_visible(False)

    fig.suptitle(
        f'All {N_SITES} Deletion Sites — Difference Maps (del − WT)\n'
        f'{CELL_NAME}  |  Each deletion = {DEL_SIZE:,} bp  |  '
        f'Colour scale ±{global_vlim:.3f}  |  Red border = Edward\'s actual insulator',
        fontsize=11, fontweight='bold'
    )

    path = 'media/deletion_scan_montage.png'
    plt.savefig(path, dpi=200, bbox_inches='tight')
    print(f'Saved: {path}')
    plt.close()
    return path


# ── Figure 4: impact ranking bar + WT triangle side by side ──────────────────
def plot_ranked_summary(records, wt):
    """
    Left: WT triangle TAD (reference).
    Right: ranked bar chart of global_abs coloured by position along genome.
    Shows clearly where the most sensitive sites are relative to the TAD structure.
    """
    sorted_by_pos  = sorted(records, key=lambda r: r['centre'])
    sorted_by_rank = sorted(records, key=lambda r: r['global_abs'], reverse=True)

    cmap_pos = plt.get_cmap('coolwarm')
    pos_norm = plt.Normalize(WIN_START, WIN_END)

    fig = plt.figure(figsize=(16, 6))
    gs  = mgs.GridSpec(1, 2, width_ratios=[1.4, 1], wspace=0.35)

    # Left: WT triangle
    ax_tri = fig.add_subplot(gs[0])
    start, end = WIN_START, WIN_END
    y_max = (end - start) / 2
    vmax  = np.percentile(wt[wt > 0], 99) if np.any(wt > 0) else 1
    cmap_obj = plt.get_cmap('Reds').copy(); cmap_obj.set_bad('white')
    ax_tri.imshow(np.ma.masked_invalid(_rotate45(wt)), cmap=cmap_obj,
                  aspect='auto', vmin=0, vmax=vmax, origin='lower',
                  extent=[start, end, 0, y_max], interpolation='nearest')
    # Mark all scan positions
    for r in sorted_by_pos:
        col = '#e74c3c' if r['is_actual'] else cmap_pos(pos_norm(r['centre']))
        ax_tri.axvline(r['centre'] + DEL_SIZE // 2, color=col,
                       lw=2 if r['is_actual'] else 0.8,
                       ls='-' if r['is_actual'] else ':',
                       alpha=1.0 if r['is_actual'] else 0.6)
    ax_tri.set_xlim(start, end); ax_tri.set_ylim(0, y_max)
    ax_tri.set_title('Wild-type TAD structure\nColoured lines = scan positions '
                     '(red = Edward\'s insulator)', fontsize=9, fontweight='bold')
    ax_tri.set_xlabel('Genomic position (bp)', fontsize=8)
    ax_tri.set_ylabel('Distance/2 (bp)', fontsize=8)
    ax_tri.ticklabel_format(style='plain', axis='both')

    # Right: ranked impact bars
    ax_bar = fig.add_subplot(gs[1])
    centres_ranked = [r['centre'] for r in sorted_by_rank]
    vals_ranked    = [r['global_abs'] for r in sorted_by_rank]
    colours_ranked = ['#e74c3c' if r['is_actual']
                      else cmap_pos(pos_norm(r['centre']))
                      for r in sorted_by_rank]
    y_pos = np.arange(len(records))
    bars  = ax_bar.barh(y_pos, vals_ranked, color=colours_ranked,
                        edgecolor='white', linewidth=0.5)
    # Thicker border for actual
    for j, r in enumerate(sorted_by_rank):
        if r['is_actual']:
            bars[j].set_edgecolor('black'); bars[j].set_linewidth(2)
    ax_bar.set_yticks(y_pos)
    ax_bar.set_yticklabels(
        [f"{r['centre']//1_000_000:.3f} Mb" for r in sorted_by_rank],
        fontsize=7.5
    )
    ax_bar.invert_yaxis()
    ax_bar.set_xlabel('Mean |Δ contact| (impact score)', fontsize=9)
    ax_bar.set_title('Impact ranking\n(top = biggest structural change)',
                     fontsize=9, fontweight='bold')
    ax_bar.axvline(0, color='black', lw=0.5)

    fig.suptitle(
        f'Deletion Scan Summary — {ED_CHROM} | {CELL_NAME}\n'
        f'All deletions = {DEL_SIZE:,} bp  |  '
        f'Red = Edward\'s actual insulator at {ed_centre // 1_000_000:.3f} Mb',
        fontsize=11, fontweight='bold'
    )

    path = 'media/deletion_scan_ranked_summary.png'
    plt.savefig(path, dpi=200, bbox_inches='tight')
    print(f'Saved: {path}')
    plt.close()
    return path


# ── HTML report ────────────────────────────────────────────────────────────────
def generate_scan_html(paths, records):
    now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
    sorted_by_rank = sorted(records, key=lambda r: r['global_abs'], reverse=True)

    table_rows = ''
    for rank, r in enumerate(sorted_by_rank, 1):
        flag  = ' ⭐ ACTUAL' if r['is_actual'] else ''
        style = 'background:#fdecea;font-weight:bold' if r['is_actual'] else ''
        table_rows += (
            f'<tr style="{style}">'
            f'<td>#{rank}</td>'
            f'<td>{r["centre"]:,}</td>'
            f'<td>{r["del_s"]:,}–{r["del_e"]:,}</td>'
            f'<td>{r["global_abs"]:.5f}{flag}</td>'
            f'<td>{r["cross_gain"]:+.5f}</td>'
            f'<td>{r["ins_weakening"]:+.5f}</td>'
            f'</tr>\n'
        )

    def img(src, cap='', w='100%'):
        return (f'<figure style="margin:0 0 16px 0">'
                f'<img src="{src}" style="width:{w};border:1px solid #ddd;border-radius:6px">'
                f'<figcaption style="font-size:.82em;color:#555;margin-top:4px">{cap}'
                f'</figcaption></figure>')

    # Pre-compute captions so no backslashes appear inside f-string expressions
    cap_summary     = ('Left: WT triangle TAD plot with all scan positions marked. '
                       'Right: sites ranked by impact score.')
    cap_sensitivity = ("Three impact metrics at each scanned position. "
                       "Red = Edward's actual insulator. Higher bars = bigger structural effect.")
    cap_gallery     = 'Triangle TAD difference views for selected deletion sites.'
    cap_montage     = f'All {N_SITES} deletion sites on a shared colour scale.'

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>Deletion Scan — Edward chr12</title>
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
<h1>Deletion Sensitivity Scan — Edward chr12</h1>
<p class="meta">
  {ED_CHROM} | {CELL_NAME} ({CELL_TYPE}) | mm10 |
  {N_SITES} deletion sites × {DEL_SIZE:,} bp each | Generated: {now}
</p>

<div class="info">
  <b>Hypothesis:</b> Deletions <em>at</em> a TAD boundary cause the largest structural
  change (TADs merge); deletions <em>within</em> a TAD body cause little change.<br>
  <b>Method:</b> Predict WT once; for each of {N_SITES} evenly-spaced positions,
  simulate a {DEL_SIZE:,}-bp deletion and compare to WT using three metrics.<br>
  <b>Edward's actual insulator:</b> {ED_CHROM}:{ED_DEL_START:,}–{ED_DEL_END:,}
  (highlighted in red throughout).
</div>

<div class="legend">
  <b>Metrics explained</b><br>
  &bull; <b>Mean |Δ contact|</b> — total reorganisation of the contact map;
    higher = deletion caused more structural change.<br>
  &bull; <b>Cross-TAD contact gain</b> — average change in contacts across the
    deletion site; positive = the two domains are merging.<br>
  &bull; <b>Insulation weakening</b> — how much the local boundary strength
    decreases; positive = boundary is lost.
</div>

<h2>Impact Ranking Table</h2>
<table>
  <tr><th>Rank</th><th>Centre (bp)</th><th>Deletion range</th>
      <th>Mean |Δ contact|</th><th>Cross-TAD gain</th><th>Ins. weakening</th></tr>
  {table_rows}
</table>

<h2>Summary: WT Structure + Impact Ranking</h2>
{img(paths['summary'], cap_summary)}

<h2>Sensitivity Profile (all 3 metrics)</h2>
{img(paths['sensitivity'], cap_sensitivity)}

<h2>Triangle TAD Gallery (selected sites)</h2>
<p style="font-size:.88em;color:#444">
  Each panel shows the <b>difference map</b> (deletion − WT) as a rotated 45°
  triangle. TADs appear as triangles; the insulator is the valley between them.
  Red border = actual insulator. Orange = highest-impact non-insulator.
  Blue = interior controls. Cyan dashed lines mark the deletion boundaries.
</p>
{img(paths['gallery'], cap_gallery)}

<h2>All {N_SITES} Sites — Difference Map Montage</h2>
<p style="font-size:.88em;color:#444">
  All scanned sites on a single shared colour scale. Compare the spatial pattern
  of contact changes across sites. Red border = Edward's actual insulator.
</p>
{img(paths['montage'], cap_montage)}

</body>
</html>
"""
    out = 'deletion_scan_report.html'
    with open(out, 'w') as fh:
        fh.write(html)
    print(f'HTML report saved: {out}')
    return out


# ── Run all plots ─────────────────────────────────────────────────────────────
if __name__ == '__main__':
    os.makedirs('media', exist_ok=True)

    paths = {
        'sensitivity': plot_sensitivity_profile(records),
        'gallery':     plot_triangle_gallery(records, wt),
        'montage':     plot_diff_montage(records, wt),
        'summary':     plot_ranked_summary(records, wt),
    }

    generate_scan_html(paths, records)

    print('\n── Results ────────────────────────────────────────────────────────')
    sorted_by_rank = sorted(records, key=lambda r: r['global_abs'], reverse=True)
    for rank, r in enumerate(sorted_by_rank, 1):
        flag = '  ← ACTUAL INSULATOR' if r['is_actual'] else ''
        print(f'  #{rank:2d}  {r["centre"]:,}  global_abs={r["global_abs"]:.5f}{flag}')

    actual = next(r for r in records if r['is_actual'])
    actual_rank = sorted_by_rank.index(actual) + 1
    print(f'\nEdward\'s insulator ranks #{actual_rank} out of {N_SITES} by global impact.')
    print('\nCheck media/ for PNGs and deletion_scan_report.html for the full report.')
