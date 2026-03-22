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
OUT_PDB     = "data/processed/pipeline_structure.pdb"
OUT_PNG     = "media/pipeline_stage3_structure.png"


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


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    load_dotenv()
    api_key = os.getenv("NVIDIA_ESM_FOLD_API_KEY") or os.getenv("NVIDIA_API_KEY")
    if not api_key:
        raise EnvironmentError("NVIDIA_ESM_FOLD_API_KEY (or NVIDIA_API_KEY) not set in .env")

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
