#!/usr/bin/env python3
"""
Stage 4 — Molecule generation (MolMIM) + docking (DiffDock) via NVIDIA NIM.

MolMIM optimises molecules from a seed SMILES for drug-likeness (QED).
DiffDock docks each generated molecule against the ESMFold PDB structure.

Correct API formats (from NVIDIA official examples):
  MolMIM  → seed SMILES in "smi" field; "algorithm", "property_name", etc.
  DiffDock → URL /biology/mit/diffdock; ligand as SDF text; protein = ATOM lines only

Estimated runtime: ~10 min (20 sequential DiffDock calls at 5–15 s each)

Outputs:
    data/processed/pipeline_molecules.json
    media/pipeline_stage4_molecules.png
    media/pipeline_stage4_molecules_2d.png   (optional, requires rdkit)
"""

import os
import json
import time

import numpy as np
import matplotlib.pyplot as plt
from dotenv import load_dotenv

# ── Constants ──────────────────────────────────────────────────────────────────
MOLMIM_URL   = "https://health.api.nvidia.com/v1/biology/nvidia/molmim/generate"
DIFFDOCK_URL = "https://health.api.nvidia.com/v1/biology/mit/diffdock"

IN_FASTA  = "data/processed/pipeline_target_sequence.fasta"
IN_PDB    = "data/processed/pipeline_structure.pdb"
OUT_JSON  = "data/processed/pipeline_molecules.json"
OUT_PNG   = "media/pipeline_stage4_molecules.png"
OUT_PNG_2D = "media/pipeline_stage4_molecules_2d.png"

NUM_MOLECULES = 20
DOCKING_POSES = 10

# Seed SMILES for MolMIM — a drug-like scaffold used as starting point for
# CMA-ES optimisation toward QED (drug-likeness). Death domains are α-helical
# bundle receptors; an aromatic/heterocyclic scaffold is a reasonable prior.
# Source: NVIDIA NIM MolMIM example (ergoline derivative, drug-like, QED ~0.7)
SEED_SMILES = "[H][C@@]12Cc3c[nH]c4cccc(C1=C[C@H](NC(=O)N(CC)CC)CN2C)c34"


# ── Helpers ────────────────────────────────────────────────────────────────────

def read_fasta(path: str = IN_FASTA) -> str:
    lines = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith(">"):
                lines.append(line)
    return "".join(lines)


def filter_pdb_atoms(pdb_string: str) -> str:
    """Return only ATOM lines from PDB (as required by DiffDock)."""
    return "\n".join(line for line in pdb_string.splitlines() if line.startswith("ATOM"))


def smiles_to_sdf(smiles: str) -> str | None:
    """Convert SMILES to 3D SDF.

    Tries rdkit first (best), then PubChem REST API (no key, free).
    Returns None on failure — caller will skip or use SMILES fallback.
    """
    # 1. rdkit
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
            AllChem.MMFFOptimizeMolecule(mol)
            return Chem.MolToMolBlock(mol)
    except Exception:
        pass

    # 2. PubChem REST API (free, 3D-optimised SDF)
    try:
        import urllib.request, urllib.parse
        encoded = urllib.parse.quote(smiles, safe="")
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{encoded}/SDF?record_type=3d"
        with urllib.request.urlopen(url, timeout=15) as r:
            sdf = r.read().decode()
        if "$$$$" in sdf:   # valid SDF marker
            return sdf
    except Exception:
        pass

    return None


def call_molmim(api_key: str, num_molecules: int = NUM_MOLECULES) -> list[dict]:
    """Generate candidate molecules via CMA-ES optimisation from seed SMILES.
    Returns list of {"smiles": ..., "score": ...} dicts.
    """
    import urllib.request

    payload = json.dumps({
        "algorithm": "CMA-ES",
        "num_molecules": num_molecules,
        "property_name": "QED",
        "minimize": False,
        "min_similarity": 0.3,
        "particles": 30,
        "iterations": 10,
        "smi": SEED_SMILES,
    }).encode("utf-8")
    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json",
        "Accept": "application/json",
    }

    print(f"Calling MolMIM (CMA-ES, QED optimisation, {num_molecules} molecules) ...")
    req = urllib.request.Request(MOLMIM_URL, data=payload, headers=headers, method="POST")
    try:
        with urllib.request.urlopen(req, timeout=120) as resp:
            body = json.loads(resp.read().decode())
    except Exception as exc:
        raise RuntimeError(f"MolMIM API call failed: {exc}")

    # Response key is "molecules" (v2) or "generated" (older docs)
    mol_list = body.get("molecules") or body.get("generated")
    if mol_list is None:
        print(f"  Unexpected response keys: {list(body.keys())}")
        raise KeyError(f"No molecules key in MolMIM response. Keys: {list(body.keys())}")

    # molecules is a JSON-encoded string — double-parse it
    if isinstance(mol_list, str):
        mol_list = json.loads(mol_list)

    print(f"  MolMIM returned {len(mol_list)} molecules.")
    # Each item: {"sample": "SMILES", "score": float}
    result = []
    for item in mol_list:
        if isinstance(item, str):
            result.append({"smiles": item, "score": 0.0})
        elif isinstance(item, dict):
            smiles = item.get("sample") or item.get("smiles", "")
            score  = float(item.get("score", 0.0))
            result.append({"smiles": smiles, "score": score})
    return result


def call_diffdock(atom_pdb: str, smiles: str, api_key: str) -> float:
    """Dock one molecule against the protein. Returns best pose_confidence (0.0 on failure).

    Ligand format: SDF preferred (via rdkit); falls back to raw SMILES string.
    Protein: ATOM-lines-only PDB string (pre-filtered).
    """
    import urllib.request

    # Prefer SDF (3D conformer); fall back to SMILES string
    sdf = smiles_to_sdf(smiles)
    if sdf:
        ligand_content = sdf
        ligand_type = "sdf"
    else:
        ligand_content = smiles
        ligand_type = "smiles"

    payload = json.dumps({
        "ligand": ligand_content,
        "ligand_file_type": ligand_type,
        "protein": atom_pdb,
        "num_poses": DOCKING_POSES,
        "time_divisions": 20,
        "steps": 18,
        "save_trajectory": False,
    }).encode("utf-8")
    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json",
        "Accept": "application/json",
    }

    import urllib.error
    try:
        req = urllib.request.Request(DIFFDOCK_URL, data=payload, headers=headers, method="POST")
        with urllib.request.urlopen(req, timeout=90) as resp:
            body = json.loads(resp.read().decode())
        # Response key confirmed: position_confidence (Array of Float32)
        for key in ("position_confidence", "pose_confidence", "confidence_scores", "confidence"):
            val = body.get(key)
            if val:
                scores = val if isinstance(val, list) else [val]
                return float(max(scores))
        return 0.0
    except urllib.error.HTTPError as exc:
        if exc.code == 429:
            print(f"    DiffDock rate-limited — waiting 10s ...")
            time.sleep(10)
        else:
            print(f"    DiffDock failed ({smiles[:30]}...): HTTP {exc.code}")
        return 0.0
    except Exception as exc:
        print(f"    DiffDock failed ({smiles[:30]}...): {exc}")
        return 0.0


def run_docking_batch(atom_pdb: str, molecules: list[dict],
                      api_key: str) -> list[dict]:
    """Run DiffDock for each molecule. Returns ranked list."""
    n = len(molecules)
    results = []
    for i, mol in enumerate(molecules):
        smiles = mol.get("smiles", "")
        molmim_score = float(mol.get("score", 0.0))
        print(f"  Docking molecule {i + 1}/{n} ...")
        confidence = call_diffdock(atom_pdb, smiles, api_key)
        results.append({
            "smiles": smiles,
            "molmim_score": molmim_score,
            "diffdock_confidence": confidence,
        })
        if i < n - 1:
            time.sleep(0.3)

    results.sort(key=lambda x: x["diffdock_confidence"], reverse=True)
    return [{"rank": i + 1, **r} for i, r in enumerate(results)]


def save_molecules_json(ranked: list[dict], out_path: str = OUT_JSON) -> None:
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(ranked, f, indent=2)
    print(f"Molecules JSON saved → {out_path}")


def plot_molecules_bar(ranked: list[dict], out_png: str = OUT_PNG) -> None:
    os.makedirs("media", exist_ok=True)
    top10  = ranked[:10]
    labels  = [f"#{m['rank']}\n{m['smiles'][:12]}..." for m in top10]
    dd_conf = [m["diffdock_confidence"] for m in top10]
    mm_score = [m["molmim_score"] for m in top10]
    x     = np.arange(len(top10))
    width = 0.38

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.bar(x - width / 2, dd_conf,  width, label="DiffDock confidence",
           color="#4527A0", alpha=0.85, edgecolor="white")
    ax.bar(x + width / 2, mm_score, width, label="MolMIM QED score",
           color="#7B1FA2", alpha=0.65, edgecolor="white")
    ax.set_xlabel("Molecule (rank · SMILES prefix)", fontsize=11)
    ax.set_ylabel("Score", fontsize=11)
    ax.set_title(
        "UNC5B Death Domain — Top-10 Drug Candidates\n"
        "DiffDock Pose Confidence + MolMIM QED Score",
        fontsize=12, fontweight="bold",
    )
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7.5, rotation=20, ha="right")
    ax.legend(fontsize=9)
    ax.set_ylim(0, 1.1)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Bar chart saved → {out_png}")


def plot_molecules_2d(ranked: list[dict], out_png: str = OUT_PNG_2D) -> None:
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw
    except ImportError:
        return
    os.makedirs("media", exist_ok=True)
    mols, legends = [], []
    for m in ranked[:10]:
        mol = Chem.MolFromSmiles(m["smiles"])
        if mol is not None:
            mols.append(mol)
            legends.append(f"#{m['rank']} conf={m['diffdock_confidence']:.2f}")
    if not mols:
        return
    img = Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(300, 200), legends=legends)
    img.save(out_png)
    print(f"2D grid saved → {out_png}")


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    load_dotenv()
    molmim_key   = os.getenv("NVIDIA_MOLMIM_API_KEY") or os.getenv("NVIDIA_API_KEY")
    diffdock_key = os.getenv("NVIDIA_DIFF_API_KEY")   or os.getenv("NVIDIA_API_KEY")
    if not molmim_key:
        raise EnvironmentError("NVIDIA_MOLMIM_API_KEY (or NVIDIA_API_KEY) not set in .env")
    if not diffdock_key:
        raise EnvironmentError("NVIDIA_DIFF_API_KEY (or NVIDIA_API_KEY) not set in .env")

    if os.path.exists(OUT_JSON):
        print(f"Molecules cache exists → {OUT_JSON} (skip API calls; delete to rerun)")
        try:
            with open(OUT_JSON) as f:
                ranked = json.load(f)
            if not isinstance(ranked, list) or not ranked:
                raise ValueError("Empty or invalid molecules JSON.")
            print(f"  Cached: {len(ranked)} molecules")
            plot_molecules_bar(ranked)
            plot_molecules_2d(ranked)
            print("Stage 4 complete.")
            return
        except Exception as exc:
            print(f"  Cache corrupt ({exc}); recomputing...")

    if not os.path.exists(IN_PDB):
        raise FileNotFoundError(f"PDB not found: {IN_PDB}. Run Stage 3 first.")

    with open(IN_PDB) as f:
        pdb_string = f.read()
    atom_pdb = filter_pdb_atoms(pdb_string)
    if not atom_pdb.strip():
        raise ValueError(f"No ATOM lines found in PDB: {IN_PDB}")
    print(f"PDB loaded: {atom_pdb.count(chr(10))} ATOM lines")

    molecules = call_molmim(molmim_key)

    print(f"Running DiffDock on {len(molecules)} molecules (estimated ~10 min)...")
    ranked = run_docking_batch(atom_pdb, molecules, diffdock_key)

    save_molecules_json(ranked)
    plot_molecules_bar(ranked)
    plot_molecules_2d(ranked)
    print("Stage 4 complete.")


if __name__ == "__main__":
    main()
