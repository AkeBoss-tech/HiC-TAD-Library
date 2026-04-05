#!/usr/bin/env python3
"""
Stage 3 — ESMFold structure prediction via NVIDIA NIM.

Folds the UNC5B death domain sequence and renders a 3D backbone trace.

Outputs:
    data/processed/pipeline_structure.pdb
    media/pipeline_stage3_structure.png
"""

import os
import json

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (registers 3D projection)
from dotenv import load_dotenv

# ── Constants ──────────────────────────────────────────────────────────────────
ESMFOLD_URL = "https://health.api.nvidia.com/v1/biology/nvidia/esmfold"
IN_FASTA    = "data/processed/pipeline_target_sequence.fasta"
IN_MULTI_FASTA = "data/processed/pipeline_target_sequences.fasta"
OUT_PDB     = "data/processed/pipeline_structure.pdb"
OUT_DIR     = "data/processed/pipeline_structures"
IN_EMBEDDING_JSON = "data/processed/pipeline_embedding_comparison.json"
OUT_MULTI_JSON = "data/processed/pipeline_structure_comparison.json"
OUT_PNG     = "media/pipeline_stage3_structure.png"
OUT_MULTI_PNG = "media/pipeline_stage3_structure_comparison.png"


# ── Helpers ────────────────────────────────────────────────────────────────────

def read_fasta(path: str = IN_FASTA) -> str:
    """Read FASTA file and return sequence string (no header lines)."""
    lines = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith(">"):
                lines.append(line)
    return "".join(lines)


def read_multi_fasta(path: str = IN_MULTI_FASTA) -> list[dict]:
    records = []
    header = None
    chunks = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append({"header": header, "sequence": "".join(chunks)})
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        records.append({"header": header, "sequence": "".join(chunks)})
    return records


def call_esmfold(sequence: str, api_key: str) -> str:
    """POST to ESMFold endpoint. Returns PDB string."""
    import urllib.request

    # ESMFold takes singular "sequence" (not an array)
    payload = json.dumps({"sequence": sequence}).encode("utf-8")
    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json",
        "Accept": "application/json",
    }

    print(f"Calling ESMFold endpoint for {len(sequence)}-aa sequence ...")
    req = urllib.request.Request(ESMFOLD_URL, data=payload, headers=headers, method="POST")
    try:
        with urllib.request.urlopen(req, timeout=120) as resp:
            body = json.loads(resp.read().decode())
    except Exception as exc:
        raise RuntimeError(f"ESMFold API call failed: {exc}")

    # Handle multiple possible response key names
    pdb_str = None
    for key in ("pdbs", "pdb_structure", "pdb", "structure"):
        val = body.get(key)
        if val:
            pdb_str = val[0] if isinstance(val, list) else val
            break

    if pdb_str is None:
        print(f"  Unexpected response keys: {list(body.keys())}")
        raise KeyError(f"No PDB key found in ESMFold response. Keys: {list(body.keys())}")

    print(f"  PDB string received: {len(pdb_str)} characters")
    return pdb_str


def parse_ca_coordinates(pdb_string: str) -> np.ndarray:
    """Extract Cα coordinates from PDB string. Returns ndarray (n_residues, 3)."""
    coords = []
    for line in pdb_string.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords.append([x, y, z])
    if not coords:
        raise ValueError("No Cα atoms found in PDB string.")
    arr = np.array(coords)
    print(f"  Parsed {len(arr)} Cα residues from PDB.")
    return arr


def plot_backbone_3d(coords: np.ndarray, out_png: str = OUT_PNG) -> None:
    """3D backbone trace colored by residue index (plasma colormap)."""
    os.makedirs("media", exist_ok=True)

    n = len(coords)
    colors = plt.cm.plasma(np.linspace(0, 1, n))

    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection="3d")

    # Backbone trace
    ax.plot(coords[:, 0], coords[:, 1], coords[:, 2],
            color="gray", lw=1.0, alpha=0.5, zorder=1)

    # Scatter colored by residue index
    sc = ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2],
                    c=np.arange(n), cmap="plasma", s=30, zorder=2)

    cbar = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.08)
    cbar.set_label("Residue index", fontsize=9)

    ax.set_title(
        "UNC5B Death Domain — ESMFold Predicted Structure\n"
        "Cα backbone · colored by residue index (N→C: purple→yellow)",
        fontsize=11, fontweight="bold",
    )
    ax.set_axis_off()

    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved → {out_png}")


def kabsch_rmsd(reference: np.ndarray, mobile: np.ndarray) -> float:
    n = min(len(reference), len(mobile))
    if n < 3:
        return float("nan")

    ref = reference[:n] - reference[:n].mean(axis=0)
    mob = mobile[:n] - mobile[:n].mean(axis=0)
    cov = mob.T @ ref
    v, _, wt = np.linalg.svd(cov)
    d = np.sign(np.linalg.det(v @ wt))
    rot = v @ np.diag([1.0, 1.0, d]) @ wt
    aligned = mob @ rot
    return float(np.sqrt(np.mean(np.sum((aligned - ref) ** 2, axis=1))))


def plot_differential_structures(results: list[dict], out_png: str = OUT_MULTI_PNG) -> None:
    os.makedirs("media", exist_ok=True)
    variants = [row for row in results if row["id"] != "wt_canonical"]
    if not variants:
        variants = results

    x = np.arange(len(variants))
    priority = [row.get("causal_priority", row["rmsd_to_wt"]) for row in variants]
    labels = [row["id"] for row in variants]
    e1_scores = [row.get("e1_score", 0.0) for row in variants]

    # Blue if active (A), Orange if collapsed/Boundary loss, Grey if B
    colors = ["#1A237E" if e > 0.1 else ("#E65100" if e == 0 else "#616161") for e in e1_scores]
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.bar(x, priority, color=colors, edgecolor="white", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("Causal Priority (RMSD × CBD)", fontsize=11)
    ax.set_xlabel("Genomic Variant / Isoform", fontsize=11)
    ax.set_title(
        "Differential Structure Scanner — WT vs Variants (including Locus E1 State)",
        fontsize=12,
        fontweight="bold",
    )
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#1A237E", label="Compartment A"),
        Patch(facecolor="#E65100", label="Boundary Loss (E1 set to 0)"),
        Patch(facecolor="#616161", label="Compartment B")
    ]
    ax.legend(handles=legend_elements, fontsize=9, loc="upper right")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved → {out_png}")


# ── Main ───────────────────────────────────────────────────────────────────────

def main(mode: str = "single"):
    load_dotenv()
    api_key = os.getenv("NVIDIA_ESM_FOLD_API_KEY") or os.getenv("NVIDIA_API_KEY")
    if not api_key:
        raise EnvironmentError("NVIDIA_ESM_FOLD_API_KEY (or NVIDIA_API_KEY) not set in .env")

    if mode == "differential":
        if os.path.exists(OUT_MULTI_JSON):
            print(f"Differential structure cache exists → {OUT_MULTI_JSON}")
            with open(OUT_MULTI_JSON) as f:
                results = json.load(f)
            plot_differential_structures(results)
            print("Stage 3 complete.")
            return

        if not os.path.exists(IN_MULTI_FASTA):
            raise FileNotFoundError(f"Multi-FASTA not found: {IN_MULTI_FASTA}. Run Stage 1 differential mode first.")

        os.makedirs(OUT_DIR, exist_ok=True)
        records = read_multi_fasta()
        coords_by_id = {}
        results = []

        for record in records:
            seq_id = record["header"].split("|", 1)[0]
            pdb_path = os.path.join(OUT_DIR, f"{seq_id}.pdb")
            if os.path.exists(pdb_path):
                with open(pdb_path) as f:
                    pdb_str = f.read()
            else:
                pdb_str = call_esmfold(record["sequence"], api_key)
                with open(pdb_path, "w") as f:
                    f.write(pdb_str)
                print(f"PDB saved → {pdb_path}")
            coords_by_id[seq_id] = parse_ca_coordinates(pdb_str)

        wt_coords = coords_by_id["wt_canonical"]
        # ── Biology-Scale Awareness Integration ──────────────────────────────────────
        embedding_data = {}
        if os.path.exists(IN_EMBEDDING_JSON):
            with open(IN_EMBEDDING_JSON) as f:
                rows = json.load(f)
                embedding_data = {r["id"]: r for r in rows}
        
        for record in records:
            seq_id = record["header"].split("|", 1)[0]
            rmsd = 0.0 if seq_id == "wt_canonical" else kabsch_rmsd(wt_coords, coords_by_id[seq_id])
            emb_info = embedding_data.get(seq_id, {})
            cbd = emb_info.get("composite_biological_distance", 1.0) # default if missing
            
            # Causal Priority Score: Multiplies structural damage by functional risk
            causal_priority = rmsd * cbd
            
            results.append(
                {
                    "id": seq_id,
                    "ca_atoms": int(len(coords_by_id[seq_id])),
                    "rmsd_to_wt": rmsd,
                    "composite_biological_distance": cbd,
                    "causal_priority": float(causal_priority),
                    "pdb_path": os.path.join(OUT_DIR, f"{seq_id}.pdb"),
                    "e1_score": emb_info.get("e1_score", 0.0)
                }
            )

        with open(OUT_MULTI_JSON, "w") as f:
            json.dump(results, f, indent=2)
        print(f"Structure comparison saved → {OUT_MULTI_JSON}")
        plot_differential_structures(results)
        print("Stage 3 complete.")
        return

    # Cache check
    if os.path.exists(OUT_PDB):
        print(f"PDB cache exists → {OUT_PDB} (skip API call; delete to rerun)")
        try:
            with open(OUT_PDB) as f:
                pdb_str = f.read()
            if not pdb_str.strip():
                raise ValueError("PDB file is empty.")
            coords = parse_ca_coordinates(pdb_str)
        except Exception as exc:
            print(f"  Cache corrupt ({exc}); recomputing...")
            pdb_str = None

        if pdb_str is not None:
            plot_backbone_3d(coords)
            print("Stage 3 complete.")
            return

    if not os.path.exists(IN_FASTA):
        raise FileNotFoundError(
            f"FASTA not found: {IN_FASTA}. Run Stage 1 first (make run-pipeline-stage1)."
        )

    sequence = read_fasta()
    print(f"Sequence loaded: {len(sequence)} aa")

    pdb_str = call_esmfold(sequence, api_key)

    os.makedirs(os.path.dirname(OUT_PDB), exist_ok=True)
    with open(OUT_PDB, "w") as f:
        f.write(pdb_str)
    print(f"PDB saved → {OUT_PDB}")

    coords = parse_ca_coordinates(pdb_str)
    plot_backbone_3d(coords)
    print("Stage 3 complete.")


if __name__ == "__main__":
    main()
