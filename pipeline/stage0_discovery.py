#!/usr/bin/env python3
"""
Stage 0 — Automated Locus Discovery.

This script scans a genomic region in a Hi-C/Micro-C file, detects the strongest 
TAD boundaries using the src/tad_boundaries.py toolkit, and identifies the 
genes most likely to be affected by boundary disruption (the 'hijacking' risk).

It generates the pipeline_variant_context.json manifest that drives the 
rest of the multi-scale pipeline (Stages 1-4).
"""

import os
import json
import sys
import cooler
import bioframe
import numpy as np
import pandas as pd
from datetime import datetime, timezone

# Add repo root to sys.path so src.* imports work
repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)

from src.tad_boundaries import multiscale_insulation, score_boundary_prominence
from src.compartments import compute_compartments

# --- Constants ---
MCOOL_PATH = "data/raw/mouse_microc.mcool"
RESOLUTION = 5000  # 5kb
REGION_OF_INTEREST = "chr10:60,000,000-62,000,000"  # 2Mb region (within Ensembl limit)
OUT_JSON = "data/processed/pipeline_variant_context.json"
ASSEMBLY = "mm10"

import requests

def get_genes_in_region(chrom, start, end, species="mus_musculus"):
    """Fetch protein-coding genes in a genomic region via Ensembl REST API."""
    # Ensembl likes chromosome '10' instead of 'chr10'
    chrom_clean = chrom.replace('chr', '')
    print(f"Fetching protein-coding genes in {chrom_clean}:{start}-{end} ({species}) via Ensembl REST...")
    # Add biotype=protein_coding to ensure we find targets for folding/docking
    url = f"https://rest.ensembl.org/overlap/region/{species}/{chrom_clean}:{start}-{end}?feature=gene;biotype=protein_coding;content-type=application/json"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        return response.json()
    except Exception as e:
        print(f"Warning: Ensembl lookup failed: {e}")
        return []

def discover_targets(mcool_path, region, resolution=RESOLUTION):
    """
    Find strongest boundaries and the genes they insulate.
    """
    print(f"Loading {mcool_path} at {resolution}bp resolution...")
    clr = cooler.Cooler(f"{mcool_path}::resolutions/{resolution}")
    
    # 1. Compute multi-scale insulation
    window_bp = 100_000
    print(f"Computing insulation scores for {region} (window={window_bp}bp)...")
    ins_table, _ = multiscale_insulation(clr, region, window_sizes_bp=[window_bp])
    
    # 2. Score boundaries by prominence
    print("Detecting and scoring boundaries...")
    scored_bounds = score_boundary_prominence(ins_table, window_bp, prominence_thresholds=(0.1, 0.3))
    
    # Filter to 'strong' boundaries
    strong_bounds = scored_bounds[scored_bounds['boundary_class'] == 'strong'].copy()
    if strong_bounds.empty:
        print("No strong boundaries found. Relaxing...")
        strong_bounds = scored_bounds[scored_bounds['boundary_class'].notna()].copy()
        
    # 3. Fetch genes via Ensembl Overlap
    chrom, s_e = region.split(':')
    start, end = map(lambda x: int(x.replace(',', '')), s_e.split('-'))
    region_genes = get_genes_in_region(chrom, start, end)
    
    # NEW: Compute A/B compartments (Cell State Context)
    print("Computing A/B compartments (Cell State)...")
    # Estimate gene density for sign orientation
    n_bins = (end - start) // resolution
    gene_density = np.zeros(n_bins)
    for g in region_genes:
        # Find which bin(s) the gene spans
        g_s = max(start, g['start'])
        g_e = min(end, g['end'])
        bin_s = (g_s - start) // resolution
        bin_e = (g_e - start) // resolution
        if 0 <= bin_s < n_bins:
            gene_density[bin_s:bin_e+1] += 1
            
    try:
        comp_df = compute_compartments(clr, region, gene_density=gene_density)
    except Exception as e:
        print(f"Warning: Compartment analysis failed: {e}")
        comp_df = pd.DataFrame()
        
    results = []
    for _, bound in strong_bounds.iterrows():
        b_mid = (bound['start'] + bound['end']) // 2
        
        # Find closest gene from REST results
        closest_gene = None
        min_dist = float('inf')
        
        for gene in region_genes:
            g_s, g_e = gene['start'], gene['end']
            g_mid = (g_s + g_e) // 2
            dist = abs(g_mid - b_mid)
            if dist < min_dist:
                min_dist = dist
                closest_gene = gene
        
        if closest_gene:
            # Get compartment for this gene (at its TSS)
            g_tss = closest_gene['start'] if closest_gene['strand'] == 1 else closest_gene['end']
            bin_idx = (g_tss - start) // resolution
            comp_status = "unknown"
            e1_score = 0.0
            if not comp_df.empty and 0 <= bin_idx < len(comp_df):
                comp_status = comp_df.iloc[bin_idx]['compartment']
                e1_score = float(comp_df.iloc[bin_idx]['E1'])
                
            results.append({
                "boundary_pos": int(b_mid),
                "prominence": float(bound['prominence']),
                "gene_symbol": closest_gene.get('external_name', closest_gene.get('id')),
                "gene_id": closest_gene['id'],
                "distance_to_boundary": int(min_dist),
                "cell_state": {
                    "compartment": comp_status,
                    "e1_score": e1_score,
                    "expression_prior": "high" if comp_status == "A" else ("low" if comp_status == "B" else "medium")
                }
            })
            
    # Sort by prominence
    results.sort(key=lambda x: x['prominence'], reverse=True)
    return results

def main(mode="single"):
    os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
    
    targets = discover_targets(MCOOL_PATH, REGION_OF_INTEREST)
    
    if not targets:
        print("Error: No suitable targets discovered.")
        return

    top_target = targets[0]
    print(f"Top target: {top_target['gene_symbol']} (State: {top_target['cell_state']['compartment']}, Prior: {top_target['cell_state']['expression_prior']})")

    # 4. Fetch Isoforms (for differential mode)
    species = "mus_musculus"
    gene_id = top_target['gene_id']
    isoforms = []
    print(f"Fetching isoforms for {gene_id}...")
    try:
        iso_url = f"https://rest.ensembl.org/lookup/id/{gene_id}?expand=1;content-type=application/json"
        iso_data = requests.get(iso_url, timeout=10).json()
        for trans in iso_data.get("Transcript", []):
            if trans.get("biotype") == "protein_coding":
                # Very simple protein sequence extraction from Ensembl
                proj_url = f"https://rest.ensembl.org/sequence/id/{trans['id']}?type=protein;content-type=application/json"
                p_resp = requests.get(proj_url, timeout=10)
                if p_resp.status_code == 200:
                    isoforms.append({
                        "transcript_id": trans['id'],
                        "display_name": trans.get('external_name', trans['id']),
                        "protein_sequence": p_resp.json()['seq']
                    })
    except Exception as e:
        print(f"Warning: Isoform fetch failed: {e}")

    # Build the context JSON
    context = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "discovery_method": "Automated Locus Scanner (Stage 0)",
        "locus": {
            "region": REGION_OF_INTEREST,
            "boundary_midpoint": top_target['boundary_pos'],
            "boundary_prominence": top_target['prominence']
        },
        "target_gene": {
            "symbol": top_target['gene_symbol'],
            "species": species,
            "gene_id": gene_id,
            "distance_to_boundary": top_target['distance_to_boundary'],
            "cell_state": top_target['cell_state']
        },
        "isoforms": isoforms,
        "interpretation": (
            f"Strong TAD boundary ({top_target['prominence']:.2f}) found near {top_target['gene_symbol']}. "
            f"Disruption poses a hijacking risk. Current Cell State: {top_target['cell_state']['compartment']} "
            f"(Expression Prior: {top_target['cell_state']['expression_prior']})."
        ),
        "pipeline_config": {
            "resolution": RESOLUTION,
            "data_source": MCOOL_PATH
        }
    }

    with open(OUT_JSON, "w") as f:
        json.dump(context, f, indent=2)
        
    print(f"Context manifest saved -> {OUT_JSON}")
    print("Ready for Stage 1 (Genomics Context).")

if __name__ == "__main__":
    main()
