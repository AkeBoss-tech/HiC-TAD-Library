#!/usr/bin/env python3
"""
Mouse insulator deletion analysis using AlphaGenome.
Regions from collaborators (mm10):
  - Jingyun: chr13:83,739,797-83,745,138
  - Edward:  chr12:27,333,532-27,336,455
Cell type: CL:0000207 (olfactory receptor cell — closest proxy to inner ear)
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from dotenv import load_dotenv

from alphagenome.data import genome, gene_annotation, transcript
from alphagenome.models import dna_client

# ── Setup ─────────────────────────────────────────────────────────────────────
load_dotenv()
API_KEY = os.getenv('ALPHA_GENOME_API_KEY')
dna_model = dna_client.create(API_KEY)
print('Model initialized.')

# ── Mouse gene annotations (mm10, GENCODE M23) ────────────────────────────────
print('Loading mouse gene annotations...')
gtf_url = (
    'https://storage.googleapis.com/alphagenome/reference/gencode/'
    'mm10/gencode.vM23.annotation.gtf.gz.feather'
)
gtf = pd.read_feather(gtf_url)
gtf_transcript = gene_annotation.filter_transcript_support_level(
    gene_annotation.filter_protein_coding(gtf), ['1', '2']
)
gtf_longest = gene_annotation.filter_to_longest_transcript(gtf_transcript)
longest_extractor = transcript.TranscriptExtractor(gtf_longest)
longest_extractor.cache_transcripts()
print('Annotations ready.')

CELL_TYPE = 'CL:0000207'   # olfactory receptor cell (best mouse proxy for inner ear)
ORGANISM  = dna_client.Organism.MUS_MUSCULUS


# ── Core analysis function ────────────────────────────────────────────────────
def run_deletion_analysis(chrom, deletion_start, deletion_end, label):
    deletion_len = deletion_end - deletion_start + 1
    print(f'\n{"="*70}')
    print(f'{label}')
    print(f'  Region : {chrom}:{deletion_start:,}-{deletion_end:,} ({deletion_len:,} bp)')
    print(f'  Cell   : {CELL_TYPE}')
    print(f'{"="*70}')

    mid    = (deletion_start + deletion_end) // 2
    interval = genome.Interval(
        chromosome=chrom,
        start=max(0, mid - 2**19),
        end=mid + 2**19,
    )
    print(f'Window: {interval}')

    transcripts = longest_extractor.extract(interval)
    print(f'Transcripts in window: {len(transcripts)}')

    # Wild-type
    print('Predicting wild-type...')
    wt = dna_model.predict_interval(
        interval=interval,
        requested_outputs={dna_client.OutputType.CONTACT_MAPS},
        ontology_terms=[CELL_TYPE],
        organism=ORGANISM,
    )
    wt_matrix = wt.contact_maps.values[:, :, 0]
    print(f'  WT shape: {wt_matrix.shape}')

    # Deletion variant
    variant = genome.Variant(
        chromosome=chrom,
        position=deletion_start - 1,          # 0-based
        reference_bases='N' * deletion_len,
        alternate_bases='N',
    )
    print('Predicting deletion...')
    del_out = dna_model.predict_variant(
        interval=interval,
        variant=variant,
        requested_outputs={dna_client.OutputType.CONTACT_MAPS},
        ontology_terms=[CELL_TYPE],
        organism=ORGANISM,
    )
    del_matrix = del_out.alternate.contact_maps.values[:, :, 0]
    diff_matrix = del_matrix - wt_matrix
    print(f'  Del shape: {del_matrix.shape}')

    _plot_main(wt_matrix, del_matrix, diff_matrix, interval,
               deletion_start, deletion_end, deletion_len, transcripts, label, chrom)
    _plot_extra(wt_matrix, del_matrix, interval,
                deletion_start, deletion_end, label, chrom)

    return wt_matrix, del_matrix, interval


# ── Plot 1: WT | Deletion | Difference ───────────────────────────────────────
def _plot_main(wt_matrix, del_matrix, diff_matrix, interval,
               deletion_start, deletion_end, deletion_len, transcripts, label, chrom):
    fig = plt.figure(figsize=(20, 10))
    gs  = fig.add_gridspec(2, 3, height_ratios=[1, 4], hspace=0.35, wspace=0.3)

    # Gene track
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
    ax_g.set_ylim(0, 1)
    ax_g.set_yticks([])
    ax_g.legend(loc='upper right', fontsize=8)
    ax_g.set_title(f'{label}  |  {chrom}:{deletion_start:,}–{deletion_end:,}  |  {CELL_TYPE}',
                   fontsize=11, fontweight='bold')
    ax_g.ticklabel_format(style='plain', axis='x')
    for sp in ['top', 'right', 'left']:
        ax_g.spines[sp].set_visible(False)

    extent = [interval.start, interval.end, interval.end, interval.start]
    vmax   = np.percentile(wt_matrix, 99)

    def _heatmap(ax, matrix, title, cmap='Reds', vmin=0, vmax_val=None):
        im = ax.imshow(matrix, cmap=cmap, aspect='auto', interpolation='nearest',
                       extent=extent, vmin=vmin,
                       vmax=vmax_val if vmax_val is not None else np.percentile(matrix, 99))
        for c in [deletion_start, deletion_end]:
            ax.axhline(c, color='cyan', lw=1.2, ls='--', alpha=0.8)
            ax.axvline(c, color='cyan', lw=1.2, ls='--', alpha=0.8)
        ax.set_title(title, fontsize=10, fontweight='bold')
        ax.ticklabel_format(style='plain')
        ax.set_xlabel('Position (bp)')
        ax.set_ylabel('Position (bp)')
        return im

    ax0 = fig.add_subplot(gs[1, 0])
    im0 = _heatmap(ax0, wt_matrix, 'Wild-type', vmax_val=vmax)
    plt.colorbar(im0, ax=ax0, label='Contact freq.', shrink=0.8)

    ax1 = fig.add_subplot(gs[1, 1])
    im1 = _heatmap(ax1, del_matrix, f'After deletion ({deletion_len:,} bp)', vmax_val=vmax)
    plt.colorbar(im1, ax=ax1, label='Contact freq.', shrink=0.8)

    ax2  = fig.add_subplot(gs[1, 2])
    vlim = np.percentile(np.abs(diff_matrix), 99)
    im2  = _heatmap(ax2, diff_matrix, 'Deletion − WT',
                    cmap='RdBu_r', vmin=-vlim, vmax_val=vlim)
    plt.colorbar(im2, ax=ax2, label='Δ Contact freq.', shrink=0.8)

    safe = f'{chrom}_{deletion_start}_{deletion_end}'
    path = f'media/mouse_deletion_{safe}.png'
    plt.savefig(path, dpi=200, bbox_inches='tight')
    print(f'Saved: {path}')
    plt.close()


# ── Plot 2: log2 ratio | virtual 4C | P(s) ───────────────────────────────────
def _plot_extra(wt_matrix, del_matrix, interval,
                deletion_start, deletion_end, label, chrom):
    n_bins   = wt_matrix.shape[0]
    bin_size = (interval.end - interval.start) / n_bins
    positions = np.linspace(interval.start, interval.end, n_bins)

    eps        = 1e-3   # larger eps avoids NaN in log2 for near-zero bins
    log2_ratio = np.log2((del_matrix + eps) / (wt_matrix + eps))
    log2_ratio = np.nan_to_num(log2_ratio, nan=0.0, posinf=0.0, neginf=0.0)

    del_center = (deletion_start + deletion_end) // 2
    del_bin    = int((del_center - interval.start) / bin_size)
    del_bin    = np.clip(del_bin, 0, n_bins - 1)
    v4c_wt     = wt_matrix[del_bin, :]
    v4c_del    = del_matrix[del_bin, :]

    max_d      = n_bins // 2
    dist_kb    = np.arange(1, max_d) * bin_size / 1e3
    ps_wt      = np.array([np.nanmean(np.diag(wt_matrix,  k)) for k in range(1, max_d)])
    ps_del     = np.array([np.nanmean(np.diag(del_matrix, k)) for k in range(1, max_d)])

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f'{label} — Extra Analysis\n{chrom}:{deletion_start:,}–{deletion_end:,}  |  {CELL_TYPE}',
                 fontsize=11, fontweight='bold')

    extent = [interval.start, interval.end, interval.end, interval.start]

    # Log2 ratio
    ax = axes[0]
    vlim = np.percentile(np.abs(log2_ratio), 99)
    im = ax.imshow(log2_ratio, cmap='RdBu_r', aspect='auto',
                   interpolation='nearest', extent=extent, vmin=-vlim, vmax=vlim)
    for c in [deletion_start, deletion_end]:
        ax.axhline(c, color='black', lw=1, ls='--', alpha=0.6)
        ax.axvline(c, color='black', lw=1, ls='--', alpha=0.6)
    plt.colorbar(im, ax=ax, label='log₂(deletion / WT)', shrink=0.8)
    ax.set_title('Log₂ Ratio Map\n(red = gained, blue = lost)', fontsize=10, fontweight='bold')
    ax.ticklabel_format(style='plain')
    ax.set_xlabel('Position (bp)')
    ax.set_ylabel('Position (bp)')

    # Virtual 4C
    ax = axes[1]
    ax.fill_between(positions, v4c_wt,  alpha=0.15, color='steelblue')
    ax.fill_between(positions, v4c_del, alpha=0.15, color='firebrick')
    ax.plot(positions, v4c_wt,  color='steelblue', lw=1.5, label='Wild-type')
    ax.plot(positions, v4c_del, color='firebrick',  lw=1.5, label='Deletion')
    ax.axvspan(deletion_start, deletion_end, alpha=0.1, color='red')
    ax.axvline(del_center, color='gray', lw=1, ls=':', alpha=0.7, label='Viewpoint')
    ax.set_title('Virtual 4C from deletion site', fontsize=10, fontweight='bold')
    ax.set_xlabel('Genomic position (bp)')
    ax.set_ylabel('Contact frequency')
    ax.legend(fontsize=8)
    ax.ticklabel_format(style='plain', axis='x')

    # P(s)
    ax = axes[2]
    ax.loglog(dist_kb, ps_wt,  color='steelblue', lw=1.5, label='Wild-type')
    ax.loglog(dist_kb, ps_del, color='firebrick',  lw=1.5, label='Deletion')
    ax.fill_between(dist_kb, np.minimum(ps_wt, ps_del),
                    np.maximum(ps_wt, ps_del), alpha=0.15, color='purple', label='Δ')
    ax.set_title('Contact frequency vs. distance P(s)', fontsize=10, fontweight='bold')
    ax.set_xlabel('Genomic distance (kb)')
    ax.set_ylabel('Mean contact frequency')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, which='both')

    plt.tight_layout()
    safe = f'{chrom}_{deletion_start}_{deletion_end}'
    path = f'media/mouse_deletion_{safe}_extra.png'
    plt.savefig(path, dpi=200, bbox_inches='tight')
    print(f'Saved: {path}')
    plt.close()


# ── Run both regions ──────────────────────────────────────────────────────────
if __name__ == '__main__':
    run_deletion_analysis(
        chrom='chr13',
        deletion_start=83_739_797,
        deletion_end=83_745_138,
        label="Jingyun's insulator deletion",
    )

    run_deletion_analysis(
        chrom='chr12',
        deletion_start=27_333_532,
        deletion_end=27_336_455,
        label="Edward's insulator deletion",
    )

    print('\nDone. Check media/ for output PNGs.')
