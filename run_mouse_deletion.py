#!/usr/bin/env python3
"""
Mouse insulator deletion analysis using AlphaGenome.

Regions (mm10) — coordinates confirmed by Jingyun and Edward:
  Jingyun: chr13:83,739,797-83,745,138  (~5.3 kb insulator)
    Full 2-TAD window: chr13:81,760,002-85,200,000 (3.4 Mb)
  Edward:  chr12:27,333,532-27,336,455   (~3 kb insulator)
    Full 2-TAD window: chr12:26,440,002-28,560,000 (2.1 Mb)

Each insulator is flanked by two TADs. Deleting the insulator (~5 kb) is
predicted to merge the two neighboring TADs into one larger domain.

Cell types compared:
  CL:0000207  — olfactory receptor cell  (primary sensory neuron, cranial placode origin)
  EFO:0004038 — mouse embryonic stem cell (suggested by PI)

Ontology notes from PI:
  UBERON:0001846 (internal ear) EXISTS in AlphaGenome — but only for CAGE output, NOT contact maps.
  CL:0000589 (cochlear inner hair cell) EXISTS for CAGE, NOT contact maps.
  GO otic placode terms (GO:1905040, GO:0030916, GO:0071599, GO:0043049) are NOT in the training data.
  UBERON:0003249 / UBERON:0003069 (otic placode) are also NOT in the training data.
  Only 8 mouse contact map tracks total — none are inner-ear specific.
"""

import os
import datetime
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from dotenv import load_dotenv

from alphagenome.data import genome, gene_annotation, transcript
from alphagenome.models import dna_client

# ── Setup ─────────────────────────────────────────────────────────────────────
load_dotenv()
dna_model = dna_client.create(os.getenv('ALPHA_GENOME_API_KEY'))
print('Model initialized.')

print('Loading mouse gene annotations (mm10, GENCODE M23)...')
gtf = pd.read_feather(
    'https://storage.googleapis.com/alphagenome/reference/gencode/'
    'mm10/gencode.vM23.annotation.gtf.gz.feather'
)
gtf_longest = gene_annotation.filter_to_longest_transcript(
    gene_annotation.filter_transcript_support_level(
        gene_annotation.filter_protein_coding(gtf), ['1', '2']
    )
)
extractor = transcript.TranscriptExtractor(gtf_longest)
extractor.cache_transcripts()
print('Annotations ready.\n')

ORGANISM = dna_client.Organism.MUS_MUSCULUS
CELL_TYPES = {
    'CL:0000207':  'Olfactory receptor cell',
    'EFO:0004038': 'Mouse embryonic stem cell',
}

REGIONS = [
    dict(chrom='chr13', deletion_start=83_739_797, deletion_end=83_745_138,
         label="Jingyun — chr13 insulator",
         tad_window_start=81_760_002, tad_window_end=85_200_000),
    dict(chrom='chr12', deletion_start=27_333_532, deletion_end=27_336_455,
         label="Edward — chr12 insulator",
         tad_window_start=26_440_002, tad_window_end=28_560_000),
]


# ── Helper: build analysis window ─────────────────────────────────────────────
# AlphaGenome only supports fixed sequence lengths: 16384, 131072, 524288, or
# 1048576 bp.  The full 2-TAD windows (3.4 Mb / 2.1 Mb) exceed the maximum, so
# we use a 1 MB (2^20) window centred on the deletion midpoint.
# The tad_window_start/tad_window_end fields in REGIONS are kept for reference
# and display in plot titles, but are not passed to the model.
def _window(chrom, deletion_start, deletion_end):
    mid = (deletion_start + deletion_end) // 2
    return genome.Interval(
        chromosome=chrom,
        start=max(0, mid - 2**19),
        end=mid + 2**19,
    )


# ── Helper: draw a contact heatmap panel ─────────────────────────────────────
def _heatmap(ax, matrix, title, extent, deletion_start, deletion_end,
             cmap='Reds', vmin=0, vmax=None, marker_color='cyan'):
    vmax = vmax if vmax is not None else np.percentile(matrix, 99)
    im = ax.imshow(matrix, cmap=cmap, aspect='auto', interpolation='nearest',
                   extent=extent, vmin=vmin, vmax=vmax)
    for c in [deletion_start, deletion_end]:
        ax.axhline(c, color=marker_color, lw=1.2, ls='--', alpha=0.8)
        ax.axvline(c, color=marker_color, lw=1.2, ls='--', alpha=0.8)
    ax.set_title(title, fontsize=9, fontweight='bold')
    ax.ticklabel_format(style='plain')
    ax.set_xlabel('Position (bp)', fontsize=8)
    ax.set_ylabel('Position (bp)', fontsize=8)
    return im


# ── Helper: rotate contact matrix 45° for triangle TAD view ──────────────────
def _rotate45(matrix):
    """
    Vectorised 45° rotation for the classic triangle Hi-C plot.

    For each upper-triangle element matrix[i, j] (j >= i):
      x_rotated = i + j        (ranges 0 … 2n-2, maps to genomic midpoint)
      y_rotated = j - i        (ranges 0 … n-1,  maps to genomic distance/2)

    Returns an (n, 2n) array where:
      - row 0  = diagonal (self-contact, distance = 0)  → plotted at y_min
      - row n-1 = maximum off-diagonal distance          → plotted at y_max
    NaN fills the lower triangle and corners with no data.
    """
    n = matrix.shape[0]
    out = np.full((n, 2 * n), np.nan)
    rows, cols = np.triu_indices(n)
    out[cols - rows, rows + cols] = matrix[rows, cols]
    return out


# ── Plot 1: WT | Deletion | Difference  (one per cell type) ──────────────────
def plot_main(wt, dm, interval, deletion_start, deletion_end,
              deletion_len, transcripts, label, chrom, cell_type, cell_name):
    diff = dm - wt
    fig = plt.figure(figsize=(20, 10))
    gs  = fig.add_gridspec(2, 3, height_ratios=[1, 4], hspace=0.35, wspace=0.3)

    ax_g = fig.add_subplot(gs[0, :])
    for tx in transcripts:
        ti = tx.transcript_interval
        ax_g.plot([ti.start, ti.end], [0.5, 0.5], color='steelblue', lw=2)
        name = tx.info.get('gene_name', '')
        if name:
            ax_g.text((ti.start + ti.end) / 2, 0.63, name,
                      ha='center', va='bottom', fontsize=7, fontweight='bold')
    ax_g.axvspan(deletion_start, deletion_end, alpha=0.25, color='red', label='Deletion')
    ax_g.set_xlim(interval.start, interval.end)
    ax_g.set_ylim(0, 1); ax_g.set_yticks([])
    ax_g.legend(loc='upper right', fontsize=8)
    ax_g.set_title(f'{label}  |  {cell_type} ({cell_name})', fontsize=11, fontweight='bold')
    ax_g.ticklabel_format(style='plain', axis='x')
    for sp in ['top', 'right', 'left']:
        ax_g.spines[sp].set_visible(False)

    extent = [interval.start, interval.end, interval.end, interval.start]
    vmax   = np.percentile(wt, 99)

    ax0 = fig.add_subplot(gs[1, 0])
    im0 = _heatmap(ax0, wt, 'Wild-type', extent, deletion_start, deletion_end, vmax=vmax)
    plt.colorbar(im0, ax=ax0, label='Contact freq.', shrink=0.8)

    ax1 = fig.add_subplot(gs[1, 1])
    im1 = _heatmap(ax1, dm, f'After deletion ({deletion_len:,} bp)',
                   extent, deletion_start, deletion_end, vmax=vmax)
    plt.colorbar(im1, ax=ax1, label='Contact freq.', shrink=0.8)

    ax2  = fig.add_subplot(gs[1, 2])
    vlim = np.percentile(np.abs(diff), 99)
    im2  = _heatmap(ax2, diff, 'Deletion − WT',
                    extent, deletion_start, deletion_end,
                    cmap='RdBu_r', vmin=-vlim, vmax=vlim, marker_color='black')
    plt.colorbar(im2, ax=ax2, label='Δ Contact freq.', shrink=0.8)

    safe = f'{chrom}_{deletion_start}_{deletion_end}_{cell_type.replace(":", "_")}'
    path = f'media/mouse_deletion_{safe}.png'
    plt.savefig(path, dpi=200, bbox_inches='tight')
    print(f'  Saved: {path}')
    plt.close()
    return path


# ── Plot 2: log2 ratio | virtual 4C | P(s) ───────────────────────────────────
def plot_extra(wt, dm, interval, deletion_start, deletion_end,
               label, chrom, cell_type, cell_name):
    n      = wt.shape[0]
    bsize  = (interval.end - interval.start) / n
    pos    = np.linspace(interval.start, interval.end, n)
    eps    = 1e-3
    log2r  = np.nan_to_num(np.log2((dm + eps) / (wt + eps)))

    cbin   = np.clip(int(((deletion_start + deletion_end) // 2 - interval.start) / bsize), 0, n - 1)
    max_d  = n // 2
    dist   = np.arange(1, max_d) * bsize / 1e3
    ps_wt  = np.array([np.nanmean(np.diag(wt, k)) for k in range(1, max_d)])
    ps_dm  = np.array([np.nanmean(np.diag(dm, k)) for k in range(1, max_d)])

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f'{label}  |  {cell_type} ({cell_name})\nExtra analysis',
                 fontsize=11, fontweight='bold')
    extent = [interval.start, interval.end, interval.end, interval.start]

    ax = axes[0]
    vlim = np.percentile(np.abs(log2r), 99)
    im = ax.imshow(log2r, cmap='RdBu_r', aspect='auto', interpolation='nearest',
                   extent=extent, vmin=-vlim, vmax=vlim)
    for c in [deletion_start, deletion_end]:
        ax.axhline(c, color='black', lw=1, ls='--', alpha=0.6)
        ax.axvline(c, color='black', lw=1, ls='--', alpha=0.6)
    plt.colorbar(im, ax=ax, label='log₂(del / WT)', shrink=0.8)
    ax.set_title('Log₂ Ratio Map\n(red=gained, blue=lost)', fontsize=10, fontweight='bold')
    ax.ticklabel_format(style='plain')
    ax.set_xlabel('Position (bp)'); ax.set_ylabel('Position (bp)')

    ax = axes[1]
    ax.fill_between(pos, wt[cbin, :], alpha=0.15, color='steelblue')
    ax.fill_between(pos, dm[cbin, :], alpha=0.15, color='firebrick')
    ax.plot(pos, wt[cbin, :], color='steelblue', lw=1.5, label='Wild-type')
    ax.plot(pos, dm[cbin, :], color='firebrick',  lw=1.5, label='Deletion')
    ax.axvspan(deletion_start, deletion_end, alpha=0.1, color='red')
    ax.axvline((deletion_start + deletion_end) // 2, color='gray', lw=1, ls=':', label='Viewpoint')
    ax.set_title('Virtual 4C from deletion site', fontsize=10, fontweight='bold')
    ax.set_xlabel('Genomic position (bp)'); ax.set_ylabel('Contact frequency')
    ax.legend(fontsize=8); ax.ticklabel_format(style='plain', axis='x')

    ax = axes[2]
    mask = (ps_wt > 0) & (ps_dm > 0)
    ax.loglog(dist[mask], ps_wt[mask], color='steelblue', lw=1.5, label='Wild-type')
    ax.loglog(dist[mask], ps_dm[mask], color='firebrick',  lw=1.5, label='Deletion')
    ax.fill_between(dist[mask], np.minimum(ps_wt[mask], ps_dm[mask]),
                    np.maximum(ps_wt[mask], ps_dm[mask]),
                    alpha=0.15, color='purple', label='Δ')
    ax.set_title('P(s): contact freq. vs. distance', fontsize=10, fontweight='bold')
    ax.set_xlabel('Genomic distance (kb)'); ax.set_ylabel('Mean contact freq.')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3, which='both')

    plt.tight_layout()
    safe = f'{chrom}_{deletion_start}_{deletion_end}_{cell_type.replace(":", "_")}'
    path = f'media/mouse_deletion_{safe}_extra.png'
    plt.savefig(path, dpi=200, bbox_inches='tight')
    print(f'  Saved: {path}')
    plt.close()
    return path


# ── Plot 3: side-by-side cell type comparison (square heatmaps) ──────────────
def plot_celltype_comparison(results, interval, deletion_start, deletion_end, label, chrom):
    """Show WT and deletion contact maps for all cell types in one figure."""
    n_types = len(results)
    fig, axes = plt.subplots(n_types, 3, figsize=(18, 6 * n_types))
    if n_types == 1:
        axes = [axes]
    fig.suptitle(f'Cell Type Comparison — {label}\n{chrom}:{deletion_start:,}–{deletion_end:,}',
                 fontsize=13, fontweight='bold')

    extent = [interval.start, interval.end, interval.end, interval.start]
    global_vmax = max(np.percentile(r['wt'], 99) for r in results.values())

    for row, (cell_type, data) in enumerate(results.items()):
        wt   = data['wt']
        dm   = data['dm']
        diff = dm - wt
        name = CELL_TYPES[cell_type]

        ax = axes[row][0]
        im = _heatmap(ax, wt, f'WT — {name}', extent,
                      deletion_start, deletion_end, vmax=global_vmax)
        plt.colorbar(im, ax=ax, shrink=0.8)
        ax.set_ylabel(f'{cell_type}\n{name}', fontsize=8, fontweight='bold')

        ax = axes[row][1]
        im = _heatmap(ax, dm, f'Deletion — {name}', extent,
                      deletion_start, deletion_end, vmax=global_vmax)
        plt.colorbar(im, ax=ax, shrink=0.8)

        ax   = axes[row][2]
        vlim = np.percentile(np.abs(diff), 99)
        im   = _heatmap(ax, diff, f'Difference — {name}', extent,
                        deletion_start, deletion_end,
                        cmap='RdBu_r', vmin=-vlim, vmax=vlim, marker_color='black')
        plt.colorbar(im, ax=ax, label='Δ Contact freq.', shrink=0.8)

    plt.tight_layout()
    safe = f'{chrom}_{deletion_start}_{deletion_end}'
    path = f'media/mouse_deletion_{safe}_celltype_comparison.png'
    plt.savefig(path, dpi=200, bbox_inches='tight')
    print(f'  Saved: {path}')
    plt.close()
    return path


# ── Plot 4: triangle TAD view (one per cell type) ────────────────────────────
def plot_triangles(wt, dm, interval, deletion_start, deletion_end,
                   label, chrom, cell_type, cell_name):
    """
    Classic rotated-45° TAD triangle view.

    The contact matrix is rotated so the diagonal (self-contacts) runs along
    the x-axis and TADs appear as upward-pointing dark triangles. The y-axis
    represents genomic distance / 2. Insulator boundaries are marked with
    vertical cyan dashed lines.

    Three rows:
      1. Wild-type
      2. After deletion
      3. Difference (del − WT)
    """
    start, end = interval.start, interval.end
    diff  = dm - wt
    vmax  = np.percentile(wt[wt > 0], 99) if np.any(wt > 0) else 1
    dvlim = np.percentile(np.abs(diff), 99)

    wt_tri   = _rotate45(wt)
    dm_tri   = _rotate45(dm)
    diff_tri = _rotate45(diff)

    # y extent: n rows × (resolution/2) per row = (end-start)/2 total
    y_max = (end - start) / 2
    extent = [start, end, 0, y_max]

    panels = [
        (wt_tri,   'Wild-type',         'Reds',   0,     vmax),
        (dm_tri,   'After deletion',    'Reds',   0,     vmax),
        (diff_tri, 'Difference (del − WT)', 'RdBu_r', -dvlim, dvlim),
    ]

    fig, axes = plt.subplots(3, 1, figsize=(14, 10),
                             gridspec_kw={'hspace': 0.55})
    fig.suptitle(
        f'TAD Triangle View — {label}\n'
        f'{cell_type} ({cell_name})\n'
        f'Deletion: {chrom}:{deletion_start:,}–{deletion_end:,}',
        fontsize=11, fontweight='bold'
    )

    for ax, (tri, title, cmap, vmin, vmax_use) in zip(axes, panels):
        cmap_obj = plt.get_cmap(cmap).copy()
        cmap_obj.set_bad('white')
        masked = np.ma.masked_invalid(tri)
        im = ax.imshow(masked, cmap=cmap_obj, aspect='auto',
                       vmin=vmin, vmax=vmax_use,
                       extent=extent, origin='lower',
                       interpolation='nearest')
        ax.axvline(deletion_start, color='cyan', lw=1.5, ls='--', alpha=0.9,
                   label='Deletion start')
        ax.axvline(deletion_end,   color='cyan', lw=1.5, ls='--', alpha=0.9,
                   label='Deletion end')
        ax.set_xlim(start, end)
        ax.set_ylim(0, y_max)
        ax.set_title(title, fontsize=9, fontweight='bold')
        ax.set_xlabel('Genomic position (bp)', fontsize=8)
        ax.set_ylabel('Genomic distance / 2 (bp)', fontsize=8)
        ax.ticklabel_format(style='plain', axis='both')
        lbl = 'Δ Contact freq.' if 'Difference' in title else 'Contact freq.'
        plt.colorbar(im, ax=ax, shrink=0.5, label=lbl)
        if 'Wild-type' in title:
            ax.legend(fontsize=7, loc='upper right')

    safe = f'{chrom}_{deletion_start}_{deletion_end}_{cell_type.replace(":", "_")}'
    path = f'media/mouse_deletion_{safe}_triangle.png'
    plt.savefig(path, dpi=200, bbox_inches='tight')
    print(f'  Saved: {path}')
    plt.close()
    return path


# ── Plot 5: triangle view comparison across all cell types ────────────────────
def plot_triangle_comparison(results, interval, deletion_start, deletion_end,
                             label, chrom):
    """
    Triangle TAD views for all cell types in one figure.
    Rows = cell types, Cols = WT | Deletion | Difference.
    """
    start, end = interval.start, interval.end
    y_max   = (end - start) / 2
    extent  = [start, end, 0, y_max]
    n_types = len(results)

    global_vmax = max(
        np.percentile(d['wt'][d['wt'] > 0], 99)
        for d in results.values() if np.any(d['wt'] > 0)
    )

    fig, axes = plt.subplots(n_types, 3, figsize=(18, 5 * n_types),
                             gridspec_kw={'hspace': 0.55, 'wspace': 0.35})
    if n_types == 1:
        axes = [axes]

    fig.suptitle(
        f'Triangle TAD Comparison — {label}\n'
        f'{chrom}:{deletion_start:,}–{deletion_end:,}',
        fontsize=13, fontweight='bold'
    )

    for row, (cell_type, data) in enumerate(results.items()):
        wt   = data['wt']
        dm   = data['dm']
        diff = dm - wt
        name = CELL_TYPES[cell_type]
        dvlim = np.percentile(np.abs(diff), 99)

        for col, (mat, title, cmap, vmin, vmax_use) in enumerate([
            (wt,   f'WT — {name}',         'Reds',   0,      global_vmax),
            (dm,   f'Deletion — {name}',   'Reds',   0,      global_vmax),
            (diff, f'Difference — {name}', 'RdBu_r', -dvlim, dvlim),
        ]):
            tri = _rotate45(mat)
            cmap_obj = plt.get_cmap(cmap).copy()
            cmap_obj.set_bad('white')
            masked = np.ma.masked_invalid(tri)
            ax = axes[row][col]
            im = ax.imshow(masked, cmap=cmap_obj, aspect='auto',
                           vmin=vmin, vmax=vmax_use,
                           extent=extent, origin='lower',
                           interpolation='nearest')
            ax.axvline(deletion_start, color='cyan', lw=1.2, ls='--', alpha=0.9)
            ax.axvline(deletion_end,   color='cyan', lw=1.2, ls='--', alpha=0.9)
            ax.set_xlim(start, end)
            ax.set_ylim(0, y_max)
            ax.set_title(title, fontsize=9, fontweight='bold')
            ax.set_xlabel('Genomic position (bp)', fontsize=8)
            if col == 0:
                ax.set_ylabel(f'{cell_type}\n{name}\n\nDistance/2 (bp)', fontsize=7,
                              fontweight='bold')
            else:
                ax.set_ylabel('Distance/2 (bp)', fontsize=8)
            ax.ticklabel_format(style='plain', axis='both')
            lbl = 'Δ Contact' if 'Difference' in title else 'Contact freq.'
            plt.colorbar(im, ax=ax, shrink=0.55, label=lbl)

    safe = f'{chrom}_{deletion_start}_{deletion_end}'
    path = f'media/mouse_deletion_{safe}_triangle_comparison.png'
    plt.savefig(path, dpi=200, bbox_inches='tight')
    print(f'  Saved: {path}')
    plt.close()
    return path


# ── HTML report generation ────────────────────────────────────────────────────
def _img(path, caption='', width='100%'):
    return (
        f'<figure style="margin:0 0 12px 0">'
        f'<img src="{path}" style="width:{width};border:1px solid #ddd;'
        f'border-radius:6px;" alt="{caption}">'
        f'<figcaption style="font-size:0.82em;color:#555;margin-top:4px">'
        f'{caption}</figcaption>'
        f'</figure>'
    )


def generate_html_report(all_results):
    """
    Write analysis_report.html — a self-contained HTML page embedding all PNGs.
    all_results is a list of dicts, one per region, containing saved file paths.
    """
    now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')

    sections = []
    for r in all_results:
        region_html = f"""
        <div class="region">
          <h2>{r['label']}</h2>
          <div class="info">
            <b>Deletion:</b> {r['chrom']}:{r['deletion_start']:,}–{r['deletion_end']:,}
            &nbsp;|&nbsp; <b>Size:</b> {r['deletion_end']-r['deletion_start']:,} bp
            &nbsp;|&nbsp; <b>Analysis window:</b> {r['window_start']:,}–{r['window_end']:,}
            (1 Mb, model limit)<br>
            <b>Full 2-TAD span:</b> {r['chrom']}:{r['tad_window_start']:,}–{r['tad_window_end']:,}
            ({(r['tad_window_end']-r['tad_window_start'])/1e6:.1f} Mb — too large for AlphaGenome API)
          </div>

          <h3>Triangle TAD View — all cell types</h3>
          <p style="font-size:0.88em;color:#444">
            The contact matrix is rotated 45°. TADs appear as <b>dark upward-pointing
            triangles</b> above the x-axis. The y-axis is genomic distance / 2. Cyan
            dashed lines mark the insulator deletion site. In the difference panel:
            <span style="color:#c0392b">red = gained contacts</span>,
            <span style="color:#2980b9">blue = lost contacts</span>.
          </p>
          {_img(r['triangle_comparison'],
                'Triangle TAD view — WT | Deletion | Difference for each cell type')}

          <h3>Square Heatmap — all cell types</h3>
          {_img(r['comparison'],
                'Square contact map — WT | Deletion | Difference for each cell type')}
        """

        for ct, ct_name in CELL_TYPES.items():
            ct_data = r['per_cell_type'].get(ct, {})
            region_html += f"""
          <h3>{ct} — {ct_name}</h3>
          <div class="grid">
            <div>
              <div class="label">Triangle TAD view (WT / Deletion / Difference)</div>
              {_img(ct_data.get('triangle',''), 'Triangle TAD view')}
            </div>
            <div>
              <div class="label">Square heatmap + gene track</div>
              {_img(ct_data.get('main',''), 'Square heatmap')}
            </div>
          </div>
          <div class="label">Log₂ ratio · Virtual 4C · P(s) curve</div>
          {_img(ct_data.get('extra',''), 'Extra analysis', width='100%')}
            """

        region_html += '</div>'
        sections.append(region_html)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Insulator Deletion Analysis — AlphaGenome</title>
  <style>
    body {{
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
      max-width: 1400px; margin: 0 auto; padding: 24px;
      background: #f0f2f5; color: #2c3e50;
    }}
    h1 {{ font-size: 1.7em; border-bottom: 3px solid #2980b9; padding-bottom: 8px; }}
    h2 {{ font-size: 1.3em; color: #1a5276; margin-top: 0; }}
    h3 {{ font-size: 1.05em; color: #2471a3; margin: 16px 0 6px; }}
    .region {{
      background: white; padding: 24px; margin: 20px 0;
      border-radius: 10px; box-shadow: 0 2px 8px rgba(0,0,0,0.10);
    }}
    .info {{
      background: #eaf4fb; padding: 10px 14px; border-radius: 6px;
      font-size: 0.88em; margin-bottom: 16px; line-height: 1.7;
    }}
    .grid {{
      display: grid; grid-template-columns: 1fr 1fr; gap: 16px;
      margin-bottom: 10px;
    }}
    .label {{ font-weight: 600; font-size: 0.85em; color: #555; margin: 8px 0 4px; }}
    .meta {{ font-size: 0.82em; color: #888; margin-bottom: 24px; }}
    .legend {{
      background: #fdfefe; border: 1px solid #d5d8dc;
      border-radius: 8px; padding: 14px; font-size: 0.85em;
      margin: 16px 0;
    }}
    .legend h4 {{ margin: 0 0 8px; }}
    .legend li {{ margin: 4px 0; }}
  </style>
</head>
<body>

<h1>Insulator Deletion Analysis — AlphaGenome Predictions</h1>
<p class="meta">
  Mouse mm10 &nbsp;|&nbsp; Generated: {now} &nbsp;|&nbsp;
  Model: AlphaGenome (MUS_MUSCULUS) &nbsp;|&nbsp;
  Window: 1 Mb centred on deletion (API max)
</p>

<div class="legend">
  <h4>How to read the triangle TAD plot</h4>
  <ul>
    <li><b>x-axis:</b> genomic position (bp)</li>
    <li><b>y-axis:</b> genomic distance / 2 — the farther from the x-axis, the more
        distant the two interacting loci are</li>
    <li><b>Dark upward triangles</b> = TADs (self-interacting domains)</li>
    <li><b>Valleys between triangles</b> = TAD boundaries / insulators</li>
    <li><b>Cyan dashed lines</b> = insulator deletion boundaries</li>
    <li>In difference panels: <b style="color:#c0392b">red = gained contacts</b>,
        <b style="color:#2980b9">blue = lost contacts</b> after deletion</li>
  </ul>
  <h4>Cell types used (no inner-ear cell type in AlphaGenome training data)</h4>
  <ul>
    <li><b>CL:0000207</b> — Olfactory receptor cell (cranial placode, closest proxy)</li>
    <li><b>EFO:0004038</b> — Mouse embryonic stem cell (suggested by PI)</li>
  </ul>
</div>

{''.join(sections)}

</body>
</html>
"""

    out = 'analysis_report.html'
    with open(out, 'w') as fh:
        fh.write(html)
    print(f'\nHTML report saved: {out}')
    return out


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    all_results = []

    for region in REGIONS:
        chrom   = region['chrom']
        ds      = region['deletion_start']
        de      = region['deletion_end']
        label   = region['label']
        dlen    = de - ds + 1
        interval = _window(chrom, ds, de)

        print(f'\n{"="*70}')
        print(f'{label}')
        print(f'  Deletion: {chrom}:{ds:,}-{de:,} ({dlen:,} bp)')
        print(f'  Window  : {interval}')
        print(f'{"="*70}')

        transcripts = extractor.extract(interval)
        print(f'  Transcripts in window: {len(transcripts)}')

        results      = {}
        per_cell_type = {}

        for cell_type, cell_name in CELL_TYPES.items():
            print(f'\n  [{cell_type}] {cell_name}')
            print('    Predicting WT...')
            wt_out = dna_model.predict_interval(
                interval=interval,
                requested_outputs={dna_client.OutputType.CONTACT_MAPS},
                ontology_terms=[cell_type],
                organism=ORGANISM,
            )
            wt = wt_out.contact_maps.values[:, :, 0]

            variant = genome.Variant(
                chromosome=chrom,
                position=ds - 1,
                reference_bases='N' * dlen,
                alternate_bases='N',
            )
            print('    Predicting deletion...')
            del_out = dna_model.predict_variant(
                interval=interval,
                variant=variant,
                requested_outputs={dna_client.OutputType.CONTACT_MAPS},
                ontology_terms=[cell_type],
                organism=ORGANISM,
            )
            dm = del_out.alternate.contact_maps.values[:, :, 0]
            results[cell_type] = {'wt': wt, 'dm': dm}

            p_main     = plot_main(wt, dm, interval, ds, de, dlen, transcripts,
                                   label, chrom, cell_type, cell_name)
            p_extra    = plot_extra(wt, dm, interval, ds, de,
                                    label, chrom, cell_type, cell_name)
            p_triangle = plot_triangles(wt, dm, interval, ds, de,
                                        label, chrom, cell_type, cell_name)
            per_cell_type[cell_type] = {
                'main': p_main, 'extra': p_extra, 'triangle': p_triangle,
            }

        p_comparison          = plot_celltype_comparison(results, interval, ds, de, label, chrom)
        p_triangle_comparison = plot_triangle_comparison(results, interval, ds, de, label, chrom)

        all_results.append({
            'label':             label,
            'chrom':             chrom,
            'deletion_start':    ds,
            'deletion_end':      de,
            'window_start':      interval.start,
            'window_end':        interval.end,
            'tad_window_start':  region['tad_window_start'],
            'tad_window_end':    region['tad_window_end'],
            'comparison':        p_comparison,
            'triangle_comparison': p_triangle_comparison,
            'per_cell_type':     per_cell_type,
        })

    generate_html_report(all_results)
    print('\nDone. Check media/ for PNGs and analysis_report.html for the full report.')
