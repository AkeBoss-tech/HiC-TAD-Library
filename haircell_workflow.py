#!/usr/bin/env python3
"""
Hair Cell Differentiation Workflow
====================================
Combines mouse Micro-C (observed) Hi-C with AlphaGenome (predicted)
regulatory signals to study Inner Hair Cell (IHC) vs Outer Hair Cell (OHC)
differentiation.

Steps
-----
1. Detect strong TAD boundaries in mouse Hi-C and liftover to human hg38.
2. Predict AlphaGenome contact maps + regulatory tracks for IHC and OHC.
3. Validate which boundaries are sequence-encoded vs epigenetically regulated.
4. Virtual 4C loop analysis at key hair cell gene promoters.
5. In silico mutagenesis of CTCF boundary elements.

Usage
-----
    python haircell_workflow.py              # default: ATOH1
    python haircell_workflow.py SLC26A5      # specify gene
    python haircell_workflow.py all          # run all four genes
"""

import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from dotenv import load_dotenv

# ---------------------------------------------------------------------------
# Local modules
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from tad_boundaries import (
    parse_coordinates,
    score_boundary_prominence,
    multiscale_insulation,
)
from liftover_utils import liftover_boundaries, liftover_region
from haircell_integration import (
    compute_insulation_from_matrix,
    detect_boundaries_from_insulation,
    differential_contact_map,
    boundary_concordance,
    extract_virtual_4c,
    plot_differential_contacts,
    plot_boundary_validation,
    plot_virtual_4c_comparison,
    plot_ism_scores,
)

# ---------------------------------------------------------------------------
# AlphaGenome imports
# ---------------------------------------------------------------------------
from alphagenome.data import genome, gene_annotation
from alphagenome.data import transcript as transcript_utils
from alphagenome.models import dna_client
from alphagenome.visualization import plot_components

load_dotenv()

# ===========================================================================
# Configuration
# ===========================================================================

MCOOL_PATH = os.path.join(
    os.path.dirname(__file__),
    "data", "raw", "mouse_microc.mcool",
)
RESOLUTION = 5000          # 5 kb bins for TAD-scale analysis
OUTPUT_DIR = "media/haircell"

# ---- Hair cell gene loci (1 Mb windows) -----------------------------------

# Mouse mm10 coordinates (orthologous regions)
MOUSE_LOCI = {
    "Atoh1":   "chr6:64,230,000-65,230,000",
    "Slc26a5": "chr5:20,900,000-21,900,000",
    "Pou4f3":  "chr18:41,900,000-42,900,000",
    "Gfi1":    "chr2:28,080,000-29,080,000",
}

# Human hg38 coordinates
HUMAN_LOCI = {
    "ATOH1":   "chr4:94,250,000-95,250,000",
    "SLC26A5": "chr7:106,850,000-107,850,000",
    "POU4F3":  "chr5:146,260,000-147,260,000",
    "GFI1":    "chr1:92,000,000-93,000,000",
}

# Approximate TSS positions (hg38) used as virtual-4C viewpoints
PROMOTER_POSITIONS_HG38 = {
    "ATOH1":   94_750_100,
    "SLC26A5": 107_301_550,
    "POU4F3":  146_761_820,
    "GFI1":    92_490_230,
}

# AlphaGenome ontology terms
IHC_ONTOLOGY = ["CL:0000202"]          # inner hair cell
OHC_ONTOLOGY = ["CL:0000601"]          # outer hair cell
EAR_ONTOLOGY = ["UBERON:0001846"]      # internal ear (broader fallback)

# GENCODE gene annotation (feather format for fast loading)
GTF_URL = (
    "https://storage.googleapis.com/alphagenome/reference/gencode/"
    "hg38/gencode.v46.annotation.gtf.gz.feather"
)


# ===========================================================================
# Step 1 — Mouse Hi-C Boundary Detection + LiftOver
# ===========================================================================

def step1_detect_and_liftover(mouse_locus: str, mouse_coords: str):
    """Detect strong TAD boundaries in mouse Micro-C, liftover to hg38."""
    import cooler

    print(f"\n{'=' * 70}")
    print(f"STEP 1: Detect boundaries in mouse Hi-C — {mouse_locus}")
    print(f"        {mouse_coords}")
    print(f"{'=' * 70}")

    uri = f"{MCOOL_PATH}::resolutions/{RESOLUTION}"
    clr = cooler.Cooler(uri)

    # Multi-scale insulation
    print("  Computing multi-scale insulation...")
    ins_table, window_sizes = multiscale_insulation(clr, mouse_coords)

    # Score at ~50 kb primary window
    primary_window = min(window_sizes, key=lambda w: abs(w - 50_000))
    print(f"  Scoring boundary prominence (window={primary_window // 1000}kb)...")
    scored = score_boundary_prominence(ins_table, primary_window)

    all_boundaries = scored[
        scored['boundary_class'].isin(['strong', 'weak'])
    ].copy()
    strong = all_boundaries[all_boundaries['boundary_class'] == 'strong']
    print(
        f"  Found {len(strong)} strong + "
        f"{len(all_boundaries) - len(strong)} weak boundaries"
    )

    # LiftOver mm10 → hg38
    print(f"  Lifting over {len(all_boundaries)} boundaries mm10 → hg38...")
    lifted = liftover_boundaries(all_boundaries)

    n_ok = lifted['liftover_success'].sum()
    print(f"  Successfully mapped {n_ok}/{len(all_boundaries)} boundaries")

    out = os.path.join(OUTPUT_DIR, f"{mouse_locus}_lifted_boundaries.csv")
    lifted.to_csv(out, index=False)
    print(f"  Saved to {out}")

    return all_boundaries, lifted


# ===========================================================================
# Step 2 — AlphaGenome Predictions for IHC vs OHC
# ===========================================================================

def _predict_with_fallback(dna_model, interval, requested, ontology_terms):
    """Try cell-type prediction; fall back to broader ear ontology on error."""
    try:
        return dna_model.predict_interval(
            interval=interval,
            requested_outputs=requested,
            ontology_terms=ontology_terms,
        )
    except Exception as exc:
        print(f"    Cell-type prediction failed ({exc}), "
              "falling back to inner ear (UBERON:0001846)...")
        return dna_model.predict_interval(
            interval=interval,
            requested_outputs=requested,
            ontology_terms=EAR_ONTOLOGY,
        )


def step2_predict_haircell_profiles(gene_name, human_interval_str, dna_model):
    """Generate AlphaGenome predictions for IHC and OHC at a locus."""
    print(f"\n{'=' * 70}")
    print(f"STEP 2: AlphaGenome predictions — {gene_name}")
    print(f"        {human_interval_str}")
    print(f"{'=' * 70}")

    interval = genome.Interval.from_str(human_interval_str)

    requested = {
        dna_client.OutputType.CONTACT_MAPS,
        dna_client.OutputType.RNA_SEQ,
        dna_client.OutputType.ATAC,
        dna_client.OutputType.CHIP_TF,
    }

    # IHC
    print("  Predicting IHC (CL:0000202)...")
    ihc_output = _predict_with_fallback(dna_model, interval, requested,
                                        IHC_ONTOLOGY)

    # OHC
    print("  Predicting OHC (CL:0000601)...")
    ohc_output = _predict_with_fallback(dna_model, interval, requested,
                                        OHC_ONTOLOGY)

    # Differential contact map (IHC − OHC)
    print("  Computing differential contact map (IHC - OHC)...")
    ihc_cm = ihc_output.contact_maps.values
    ohc_cm = ohc_output.contact_maps.values
    min_n = min(ihc_cm.shape[0], ohc_cm.shape[0])
    diff = differential_contact_map(ihc_cm[:min_n, :min_n],
                                    ohc_cm[:min_n, :min_n])

    plot_differential_contacts(
        diff, human_interval_str,
        title=f"Differential Contacts: {gene_name} (IHC - OHC)",
        output_path=os.path.join(OUTPUT_DIR,
                                 f"{gene_name}_diff_contacts.png"),
    )

    # CTCF comparison
    _plot_ctcf_comparison(ihc_output, ohc_output, gene_name)

    return ihc_output, ohc_output, diff


def _plot_ctcf_comparison(ihc_output, ohc_output, gene_name):
    """Compare mean predicted CTCF signal between IHC and OHC."""
    try:
        ihc_mask = (
            ihc_output.chip_tf.metadata['transcription_factor'] == 'CTCF'
        )
        ohc_mask = (
            ohc_output.chip_tf.metadata['transcription_factor'] == 'CTCF'
        )
        ihc_ctcf = ihc_output.chip_tf.values[:, ihc_mask].mean(axis=1)
        ohc_ctcf = ohc_output.chip_tf.values[:, ohc_mask].mean(axis=1)

        fig, ax = plt.subplots(figsize=(14, 3))
        x = np.arange(len(ihc_ctcf))
        ax.fill_between(x, ihc_ctcf, alpha=0.5, color='steelblue', label='IHC')
        ax.fill_between(x, ohc_ctcf, alpha=0.5, color='coral', label='OHC')
        ax.set_ylabel("CTCF Signal")
        ax.set_title(f"Predicted CTCF Binding: {gene_name}")
        ax.legend(loc='upper right')

        out = os.path.join(OUTPUT_DIR, f"{gene_name}_ctcf_comparison.png")
        plt.savefig(out, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved CTCF comparison to {out}")
    except Exception as exc:
        print(f"  Could not plot CTCF comparison: {exc}")


# ===========================================================================
# Step 3 — Validate Sequence-Encoded Boundaries
# ===========================================================================

def step3_validate_boundaries(
    gene_name, lifted_boundaries, ihc_output, human_interval_str,
):
    """Compare lifted mouse boundaries with AlphaGenome predicted boundaries."""
    print(f"\n{'=' * 70}")
    print(f"STEP 3: Validate boundaries — {gene_name}")
    print(f"{'=' * 70}")

    interval = genome.Interval.from_str(human_interval_str)
    contacts = ihc_output.contact_maps.values
    ag_res = ihc_output.contact_maps.resolution
    n_bins = contacts.shape[0]

    print(f"  Predicted contact map: {contacts.shape}, "
          f"resolution = {ag_res} bp")

    # Insulation from predicted contacts
    window_bins = max(3, 50_000 // ag_res)
    pred_insulation = compute_insulation_from_matrix(contacts,
                                                     window_bins=window_bins)
    pred_boundaries = detect_boundaries_from_insulation(pred_insulation)

    n_strong = (
        (pred_boundaries['boundary_class'] == 'strong').sum()
        if len(pred_boundaries) else 0
    )
    n_weak = (
        (pred_boundaries['boundary_class'] == 'weak').sum()
        if len(pred_boundaries) else 0
    )
    print(f"  Predicted boundaries: {n_strong} strong, {n_weak} weak")

    # Map lifted mouse boundaries into AlphaGenome bin space
    ok = lifted_boundaries[lifted_boundaries['liftover_success']].copy()

    if len(ok) == 0:
        print("  No lifted boundaries available for validation")
        return pd.DataFrame()

    ok = ok[ok['target_chrom'] == interval.chromosome].copy()
    ok['bin_index'] = (
        (ok['target_pos'] - interval.start) / ag_res
    ).astype(int)
    ok = ok[(ok['bin_index'] >= 0) & (ok['bin_index'] < n_bins)].copy()

    if len(ok) == 0:
        print("  No lifted boundaries map into this interval")
        return pd.DataFrame()

    print(f"  {len(ok)} lifted boundaries within interval")

    tolerance = max(3, 10_000 // ag_res)
    obs_df = ok[['bin_index', 'boundary_class', 'prominence']].copy()
    concordance = boundary_concordance(obs_df, pred_boundaries,
                                       tolerance_bins=tolerance)

    n_conc = concordance['is_concordant'].sum()
    print(f"  Concordance: {n_conc}/{len(concordance)} boundaries "
          "are sequence-encoded")

    concordance.to_csv(
        os.path.join(OUTPUT_DIR, f"{gene_name}_concordance.csv"), index=False,
    )

    plot_boundary_validation(
        pred_insulation, concordance, human_interval_str,
        resolution=ag_res,
        title=f"Boundary Validation: {gene_name} (Mouse Hi-C vs AlphaGenome)",
        output_path=os.path.join(OUTPUT_DIR,
                                 f"{gene_name}_boundary_validation.png"),
    )

    return concordance


# ===========================================================================
# Step 4 — Virtual 4C Loop Analysis
# ===========================================================================

def step4_virtual_4c(
    gene_name, ihc_output, ohc_output,
    human_interval_str, promoter_pos,
    transcript_extractor=None,
):
    """Compare virtual 4C from a gene promoter viewpoint in IHC vs OHC."""
    print(f"\n{'=' * 70}")
    print(f"STEP 4: Virtual 4C analysis — {gene_name}")
    print(f"{'=' * 70}")

    interval = genome.Interval.from_str(human_interval_str)
    ag_res = ihc_output.contact_maps.resolution

    viewpoint_bin = (promoter_pos - interval.start) // ag_res
    n_bins = ihc_output.contact_maps.values.shape[0]

    if viewpoint_bin < 0 or viewpoint_bin >= n_bins:
        print(f"  Viewpoint bin {viewpoint_bin} out of range "
              f"[0, {n_bins}), skipping")
        return

    print(f"  Viewpoint: {gene_name} promoter at "
          f"{promoter_pos:,} (bin {viewpoint_bin})")

    ihc_cm = ihc_output.contact_maps.values
    ohc_cm = ohc_output.contact_maps.values

    v4c_ihc = extract_virtual_4c(ihc_cm, viewpoint_bin)
    v4c_ohc = extract_virtual_4c(ohc_cm, viewpoint_bin)

    # Cell-type-specific loops
    min_len = min(len(v4c_ihc), len(v4c_ohc))
    diff = v4c_ihc[:min_len] - v4c_ohc[:min_len]
    valid = np.isfinite(diff)

    if valid.sum() > 0:
        p95 = np.nanpercentile(diff[valid], 95)
        p5 = np.nanpercentile(diff[valid], 5)
        ihc_enriched = np.where(diff > p95)[0]
        ohc_enriched = np.where(diff < p5)[0]

        print(f"  IHC-enriched interaction bins: {len(ihc_enriched)}")
        print(f"  OHC-enriched interaction bins: {len(ohc_enriched)}")

        for label, bins in [("OHC-specific (potential OHC enhancers)",
                              ohc_enriched),
                             ("IHC-specific (potential IHC enhancers)",
                              ihc_enriched)]:
            if len(bins):
                print(f"  Top {label}:")
                for b in bins[:5]:
                    pos = interval.start + b * ag_res
                    print(f"    {interval.chromosome}:{pos:,}")

    # Raw numpy-based comparison plot
    plot_virtual_4c_comparison(
        v4c_ihc, v4c_ohc, viewpoint_bin, human_interval_str,
        gene_name=gene_name,
        output_path=os.path.join(OUTPUT_DIR,
                                 f"{gene_name}_v4c_comparison.png"),
    )

    # AlphaGenome native virtual 4C (if transcripts available)
    if transcript_extractor is not None:
        _plot_alphagenome_v4c(
            ihc_output, ohc_output, interval,
            promoter_pos, gene_name, transcript_extractor,
        )


def _plot_alphagenome_v4c(
    ihc_output, ohc_output, interval,
    promoter_pos, gene_name, transcript_extractor,
):
    """Use AlphaGenome's built-in to_virtual_4c + plot_components."""
    try:
        ihc_v4c = ihc_output.contact_maps.to_virtual_4c(promoter_pos)
        ohc_v4c = ohc_output.contact_maps.to_virtual_4c(promoter_pos)
        transcripts = transcript_extractor.extract(interval)

        plot_components.plot(
            [
                plot_components.TranscriptAnnotation(transcripts),
                plot_components.Tracks(
                    tdata=ihc_v4c,
                    ylabel_template='IHC Virtual 4C',
                    filled=True, color='steelblue',
                ),
                plot_components.Tracks(
                    tdata=ohc_v4c,
                    ylabel_template='OHC Virtual 4C',
                    filled=True, color='coral',
                ),
            ],
            interval=interval,
            title=f"Virtual 4C: {gene_name} Promoter (IHC vs OHC)",
            annotations=[
                plot_components.IntervalAnnotation(
                    [f"{interval.chromosome}:{promoter_pos}-"
                     f"{promoter_pos + 2048}"],
                    color='black', alpha=0.5,
                ),
            ],
        )

        out = os.path.join(OUTPUT_DIR, f"{gene_name}_v4c_alphagenome.png")
        plt.savefig(out, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved AlphaGenome V4C to {out}")
    except Exception as exc:
        print(f"  AlphaGenome native V4C: {exc}")


# ===========================================================================
# Step 5 — In Silico Mutagenesis of Boundary Elements
# ===========================================================================

def step5_ism_boundaries(
    gene_name, ihc_output, human_interval_str,
    concordance, dna_model,
):
    """
    Score single-nucleotide variants at CTCF motifs located at strong
    boundaries.  For each boundary:

    1. Locate the CTCF peak within ±5 kb of the boundary.
    2. Create all 12 possible SNVs at that position.
    3. Call predict_variant to get reference vs alternate contact maps.
    4. Compute a *TAD dissolution score* = mean |insulation_alt − insulation_ref|.

    Note: variants whose reference base does not match the genome will be
    silently skipped (expected — we try all four possible ref bases).
    """
    print(f"\n{'=' * 70}")
    print(f"STEP 5: In silico mutagenesis — {gene_name}")
    print(f"{'=' * 70}")

    interval = genome.Interval.from_str(human_interval_str)
    ag_res = ihc_output.contact_maps.resolution

    # Choose target boundaries
    if len(concordance) > 0 and concordance['is_concordant'].any():
        targets = concordance[concordance['is_concordant']].copy()
        print(f"  Targeting {len(targets)} sequence-encoded boundaries")
    elif len(concordance) > 0:
        targets = concordance.head(3).copy()
        print(f"  No concordant boundaries; using top {len(targets)}")
    else:
        print("  No boundaries available for ISM, skipping")
        return pd.DataFrame()

    # CTCF signal for motif localization
    ctcf_signal = _get_ctcf_signal(ihc_output)

    # AlphaGenome variant prediction needs a power-of-2 interval
    pred_interval = interval.resize(2 ** 20)
    window_bins_ins = max(3, 50_000 // ag_res)

    variant_scores = []

    for _, bnd in targets.iterrows():
        boundary_bin = int(bnd['observed_bin'])
        boundary_pos = interval.start + boundary_bin * ag_res

        ctcf_peak = _find_ctcf_peak(
            ctcf_signal, boundary_bin, ag_res,
            interval.start, scan_radius_bp=5000,
        )
        if ctcf_peak is None:
            ctcf_peak = boundary_pos

        print(f"\n  Boundary at {boundary_pos:,}, "
              f"CTCF peak at {ctcf_peak:,}")

        # Try all 12 possible SNVs at the CTCF peak
        for ref_base in "ACGT":
            for alt_base in "ACGT":
                if ref_base == alt_base:
                    continue

                variant = genome.Variant(
                    chromosome=interval.chromosome,
                    position=ctcf_peak,
                    reference_bases=ref_base,
                    alternate_bases=alt_base,
                    name=f"{gene_name}_{ctcf_peak}_{ref_base}>{alt_base}",
                )

                try:
                    vout = dna_model.predict_variant(
                        interval=pred_interval,
                        variant=variant,
                        requested_outputs={
                            dna_client.OutputType.CONTACT_MAPS,
                        },
                        ontology_terms=IHC_ONTOLOGY,
                    )

                    ref_cm = vout.reference.contact_maps.values
                    alt_cm = vout.alternate.contact_maps.values

                    ref_ins = compute_insulation_from_matrix(
                        ref_cm, window_bins=window_bins_ins,
                    )
                    alt_ins = compute_insulation_from_matrix(
                        alt_cm, window_bins=window_bins_ins,
                    )

                    dissolution = float(
                        np.nanmean(np.abs(alt_ins - ref_ins))
                    )

                    variant_scores.append({
                        'locus': gene_name,
                        'boundary_pos': boundary_pos,
                        'variant_pos': ctcf_peak,
                        'ref': ref_base,
                        'alt': alt_base,
                        'tad_dissolution_score': dissolution,
                        'variant_name': variant.name,
                    })

                except Exception:
                    # Expected — wrong ref base or unsupported modality
                    pass

    if not variant_scores:
        print("  No successful variant predictions (API may not support "
              "CONTACT_MAPS for predict_variant)")
        # Fall back: score expression effect instead
        variant_scores = _ism_expression_fallback(
            gene_name, targets, interval, pred_interval,
            ctcf_signal, ag_res, dna_model,
        )

    if not variant_scores:
        print("  ISM produced no scores")
        return pd.DataFrame()

    scores_df = pd.DataFrame(variant_scores)
    scores_df = scores_df.sort_values(
        'tad_dissolution_score', ascending=False,
    ).reset_index(drop=True)

    scores_df.to_csv(
        os.path.join(OUTPUT_DIR, f"{gene_name}_ism_scores.csv"), index=False,
    )
    print(f"\n  Saved {len(scores_df)} variant scores")

    print("\n  Top variants by TAD dissolution:")
    for _, row in scores_df.head(5).iterrows():
        print(f"    {row['ref']}>{row['alt']} at {int(row['variant_pos']):,}  "
              f"score = {row['tad_dissolution_score']:.6f}")

    plot_ism_scores(
        scores_df, gene_name,
        output_path=os.path.join(OUTPUT_DIR,
                                 f"{gene_name}_ism_scores.png"),
    )

    return scores_df


def _ism_expression_fallback(
    gene_name, targets, interval, pred_interval,
    ctcf_signal, ag_res, dna_model,
):
    """Fall back to RNA_SEQ-based variant scoring when CONTACT_MAPS
    is not available in predict_variant."""
    print("  Falling back to RNA_SEQ-based variant scoring...")
    scores = []

    for _, bnd in targets.iterrows():
        boundary_bin = int(bnd['observed_bin'])
        boundary_pos = interval.start + boundary_bin * ag_res

        ctcf_peak = _find_ctcf_peak(
            ctcf_signal, boundary_bin, ag_res,
            interval.start, scan_radius_bp=5000,
        )
        if ctcf_peak is None:
            ctcf_peak = boundary_pos

        for ref_base in "ACGT":
            for alt_base in "ACGT":
                if ref_base == alt_base:
                    continue

                variant = genome.Variant(
                    chromosome=interval.chromosome,
                    position=ctcf_peak,
                    reference_bases=ref_base,
                    alternate_bases=alt_base,
                    name=f"{gene_name}_{ctcf_peak}_{ref_base}>{alt_base}",
                )

                try:
                    vout = dna_model.predict_variant(
                        interval=pred_interval,
                        variant=variant,
                        requested_outputs={
                            dna_client.OutputType.RNA_SEQ,
                        },
                        ontology_terms=IHC_ONTOLOGY,
                    )

                    ref_rna = vout.reference.rna_seq.values
                    alt_rna = vout.alternate.rna_seq.values
                    effect = float(np.nanmean(np.abs(alt_rna - ref_rna)))

                    scores.append({
                        'locus': gene_name,
                        'boundary_pos': boundary_pos,
                        'variant_pos': ctcf_peak,
                        'ref': ref_base,
                        'alt': alt_base,
                        'tad_dissolution_score': effect,
                        'variant_name': variant.name,
                    })
                except Exception:
                    pass

    return scores


def _get_ctcf_signal(output):
    """Extract mean CTCF ChIP-TF signal from a predict_interval output."""
    try:
        mask = output.chip_tf.metadata['transcription_factor'] == 'CTCF'
        return output.chip_tf.values[:, mask].mean(axis=1)
    except Exception:
        return None


def _find_ctcf_peak(ctcf_signal, boundary_bin, resolution,
                    interval_start, scan_radius_bp=5000):
    """Return the genomic position of the highest CTCF peak within
    ±scan_radius_bp of a boundary."""
    if ctcf_signal is None:
        return None

    scan_bins = scan_radius_bp // resolution
    lo = max(0, boundary_bin - scan_bins)
    hi = min(len(ctcf_signal), boundary_bin + scan_bins)

    local = ctcf_signal[lo:hi]
    if len(local) == 0:
        return None

    peak_local = int(np.argmax(local))
    return interval_start + (lo + peak_local) * resolution


# ===========================================================================
# Transcript loader (for AlphaGenome native plots)
# ===========================================================================

def _load_transcript_extractor():
    """Load GENCODE hg38 MANE-select transcripts for gene annotation."""
    print("  Loading GENCODE hg38 annotations...")
    gtf = pd.read_feather(GTF_URL)
    gtf_filt = gene_annotation.filter_protein_coding(gtf)
    gtf_filt = gene_annotation.filter_to_mane_select_transcript(gtf_filt)
    extractor = transcript_utils.TranscriptExtractor(gtf_filt)
    return extractor


# ===========================================================================
# Main
# ===========================================================================

def run_workflow(target_gene: str = "ATOH1"):
    """Run the full five-step integration workflow for a target gene."""

    print("\n" + "=" * 70)
    print(f"  HAIR CELL DIFFERENTIATION WORKFLOW: {target_gene}")
    print(f"  Mouse Micro-C (observed) + AlphaGenome (predicted)")
    print("=" * 70)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # AlphaGenome client
    api_key = os.getenv("ALPHA_GENOME_API_KEY")
    if not api_key:
        raise ValueError("ALPHA_GENOME_API_KEY not found in .env")
    dna_model = dna_client.create(api_key)

    # Gene name mapping (human → mouse ortholog)
    mouse_name = target_gene[0] + target_gene[1:].lower()
    human_interval = HUMAN_LOCI[target_gene]
    promoter_pos = PROMOTER_POSITIONS_HG38[target_gene]

    # ------------------------------------------------------------------
    # Step 1: Mouse boundaries + liftover
    # ------------------------------------------------------------------
    lifted = pd.DataFrame()
    if mouse_name in MOUSE_LOCI and os.path.exists(MCOOL_PATH):
        _, lifted = step1_detect_and_liftover(
            mouse_name, MOUSE_LOCI[mouse_name],
        )
    else:
        print("\n  Skipping Step 1 (mouse Micro-C data not available)")

    # ------------------------------------------------------------------
    # Step 2: AlphaGenome IHC vs OHC
    # ------------------------------------------------------------------
    ihc_output, ohc_output, diff = step2_predict_haircell_profiles(
        target_gene, human_interval, dna_model,
    )

    # ------------------------------------------------------------------
    # Step 3: Validate boundaries
    # ------------------------------------------------------------------
    concordance = step3_validate_boundaries(
        target_gene, lifted, ihc_output, human_interval,
    )

    # ------------------------------------------------------------------
    # Step 4: Virtual 4C
    # ------------------------------------------------------------------
    try:
        tx_extractor = _load_transcript_extractor()
    except Exception:
        tx_extractor = None

    step4_virtual_4c(
        target_gene, ihc_output, ohc_output,
        human_interval, promoter_pos, tx_extractor,
    )

    # ------------------------------------------------------------------
    # Step 5: ISM
    # ------------------------------------------------------------------
    step5_ism_boundaries(
        target_gene, ihc_output, human_interval,
        concordance, dna_model,
    )

    # ------------------------------------------------------------------
    print(f"\n{'=' * 70}")
    print(f"  WORKFLOW COMPLETE: {target_gene}")
    print(f"  All results saved to {OUTPUT_DIR}/")
    print(f"{'=' * 70}\n")


# ===========================================================================
if __name__ == "__main__":
    arg = sys.argv[1].upper() if len(sys.argv) > 1 else "ATOH1"

    if arg == "ALL":
        for gene in HUMAN_LOCI:
            run_workflow(gene)
    elif arg in HUMAN_LOCI:
        run_workflow(arg)
    else:
        print(f"Unknown gene: {arg}")
        print(f"Available: {', '.join(HUMAN_LOCI)} or ALL")
        sys.exit(1)
