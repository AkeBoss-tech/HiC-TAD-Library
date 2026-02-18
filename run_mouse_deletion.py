#!/usr/bin/env python3
"""
Mouse insulator deletion analysis using AlphaGenome.

Regions (mm10):
  Jingyun: chr13:83,739,797-83,745,138  (Mef2c locus — NOTE: PI flagged coordinates may be off)
  Edward:  chr12:27,333,532-27,336,455

Cell types compared:
  CL:0000207  — olfactory receptor cell  (primary sensory neuron, cranial placode origin)
  EFO:0004038 — mouse embryonic stem cell (suggested by PI)

Ontology notes from PI:
  UBERON:0001846 (internal ear) EXISTS in AlphaGenome — but only for CAGE output, NOT contact maps.
  GO otic placode terms (GO:1905040, GO:0030916, GO:0071599, GO:0043049) are NOT in the training data.
  UBERON:0003249 / UBERON:0003069 (otic placode) are also NOT in the training data.
  Only 8 mouse contact map tracks total — none are inner-ear specific.
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
         label="Jingyun — chr13 (Mef2c locus)"),
    dict(chrom='chr12', deletion_start=27_333_532, deletion_end=27_336_455,
         label="Edward — chr12"),
]


# ── Helper: build 1 MB window ─────────────────────────────────────────────────
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


# ── Plot 3: side-by-side cell type comparison (WT only) ──────────────────────
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


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
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

        results = {}
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

            plot_main(wt, dm, interval, ds, de, dlen, transcripts,
                      label, chrom, cell_type, cell_name)
            plot_extra(wt, dm, interval, ds, de, label, chrom, cell_type, cell_name)

        plot_celltype_comparison(results, interval, ds, de, label, chrom)

    print('\nDone. Check media/ for output PNGs.')
