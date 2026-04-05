#!/usr/bin/env python3
"""
Stage 2 — ESM-2 protein embeddings via NVIDIA NIM.

Calls the ESM2-650M embeddings endpoint with the UNC5B death domain sequence.
Response is a binary NPZ file (not JSON) — loaded with numpy and saved as .npy.

Outputs:
    data/processed/pipeline_esm2_embeddings.npy
    media/pipeline_stage2_embeddings.png
"""

import io
import os
import json
import time

import numpy as np
import matplotlib.pyplot as plt
from dotenv import load_dotenv

# ── Constants ──────────────────────────────────────────────────────────────────
ESM2_URL = "https://health.api.nvidia.com/v1/biology/meta/esm2-650m"
IN_FASTA = "data/processed/pipeline_target_sequence.fasta"
IN_MULTI_FASTA = "data/processed/pipeline_target_sequences.fasta"
IN_VARIANT_CONTEXT = "data/processed/pipeline_variant_context.json"
OUT_NPY  = "data/processed/pipeline_esm2_embeddings.npy"
OUT_MULTI_NPZ = "data/processed/pipeline_esm2_embeddings_multi.npz"
OUT_MULTI_JSON = "data/processed/pipeline_embedding_comparison.json"
OUT_PNG  = "media/pipeline_stage2_embeddings.png"
OUT_MULTI_PNG = "media/pipeline_stage2_embeddings_differential.png"


# ── Helpers ────────────────────────────────────────────────────────────────────

def read_fasta(path: str = IN_FASTA) -> str:
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


def call_esm2_embeddings(sequence: str, api_key: str) -> np.ndarray:
    """POST to ESM2-650M. Returns ndarray extracted from binary NPZ response."""
    import urllib.request

    # format="npz" → binary NPZ response body (not JSON)
    payload = json.dumps({"sequences": [sequence], "format": "npz"}).encode("utf-8")
    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json",
        "Accept": "application/json",
    }

    for attempt in (1, 2):
        try:
            req = urllib.request.Request(ESM2_URL, data=payload, headers=headers, method="POST")
            with urllib.request.urlopen(req, timeout=60) as resp:
                content_type = resp.headers.get("Content-Type", "")
                body_bytes = resp.read()
        except Exception as exc:
            if attempt == 1:
                print(f"  ESM2 call failed (attempt {attempt}): {exc}. Retrying in 5s...")
                time.sleep(5)
                continue
            raise

        print(f"  Response Content-Type: {content_type}  ({len(body_bytes)} bytes)")

        # Binary NPZ / ZIP response
        if "json" not in content_type:
            try:
                # Handle ZIP wrapping
                if content_type == "application/zip" or body_bytes[:2] == b"PK":
                    import zipfile
                    with zipfile.ZipFile(io.BytesIO(body_bytes)) as zf:
                        npz_bytes = zf.read(zf.namelist()[0])
                else:
                    npz_bytes = body_bytes
                loaded = np.load(io.BytesIO(npz_bytes), allow_pickle=False)
                # np.load returns ndarray for .npy, NpzFile (dict-like) for .npz
                if isinstance(loaded, np.ndarray):
                    emb = loaded
                    print(f"  NPY array → shape {emb.shape}")
                else:
                    for key in ("embeddings", "mean_representations", "representations", "output"):
                        if key in loaded:
                            emb = loaded[key]
                            break
                    else:
                        key = list(loaded.keys())[0]
                        emb = loaded[key]
                    print(f"  NPZ key '{key}' → shape {emb.shape}")
                return emb.astype(float)
            except Exception as exc:
                raise RuntimeError(f"Failed to parse NPZ response: {exc}")

        # Fallback: JSON response
        body = json.loads(body_bytes.decode())
        if "embeddings" not in body:
            print(f"  Unexpected JSON keys: {list(body.keys())}")
            raise KeyError(f"'embeddings' not found. Keys: {list(body.keys())}")
        raw = body["embeddings"]
        if isinstance(raw, list) and raw and isinstance(raw[0], dict):
            emb = np.array(raw[0]["embedding"], dtype=float)
        else:
            emb = np.array(raw[0] if isinstance(raw[0], list) else raw, dtype=float)
            if emb.ndim == 2 and emb.shape[0] > 2:
                emb = emb[1:-1]
        print(f"  ESM2 embedding shape: {emb.shape}")
        return emb

    raise RuntimeError("ESM2 call failed after 2 attempts.")


def plot_embeddings(embeddings: np.ndarray, out_png: str = OUT_PNG) -> None:
    """Bar chart of top-20 embedding dimensions by |value|."""
    os.makedirs("media", exist_ok=True)

    # Flatten to 1D: mean-pool if per-residue (L, D), else use as-is (D,)
    mean_emb = embeddings.flatten() if embeddings.ndim == 1 else embeddings.mean(axis=0)
    n_dims  = len(mean_emb)
    n_show  = min(20, n_dims)
    top_idx = np.argsort(np.abs(mean_emb))[-n_show:][::-1]
    values  = mean_emb[top_idx]
    colors  = ["steelblue" if v >= 0 else "salmon" for v in values]

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(range(n_show), values, color=colors, edgecolor="white", linewidth=0.5)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(range(n_show))
    ax.set_xticklabels([f"dim {i}" for i in top_idx], rotation=45, ha="right", fontsize=8)
    ax.set_xlabel("ESM-2 embedding dimension", fontsize=11)
    ax.set_ylabel("Value", fontsize=11)
    ax.set_title(
        f"UNC5B Death Domain — ESM-2 650M Embeddings ({n_dims}-dim)\n"
        f"Top-{n_show} dimensions by |value|",
        fontsize=12, fontweight="bold",
    )
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor="steelblue", label="Positive"),
                       Patch(facecolor="salmon",    label="Negative")],
              fontsize=9, loc="upper right")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved → {out_png}")


def mean_pool_embedding(embedding: np.ndarray) -> np.ndarray:
    return embedding.astype(float) if embedding.ndim == 1 else embedding.mean(axis=0).astype(float)


def cosine_distance(a: np.ndarray, b: np.ndarray) -> float:
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denom == 0:
        return 0.0
    similarity = float(np.dot(a, b) / denom)
    return 1.0 - similarity


def plot_differential_embeddings(results: list[dict], out_png: str = OUT_MULTI_PNG) -> None:
    os.makedirs("media", exist_ok=True)
    variants = [row for row in results if row["id"] != "wt_canonical"]
    if not variants:
        variants = results

    x = np.arange(len(variants))
    distances = [row.get("composite_biological_distance", row["cosine_distance_to_wt"]) for row in variants]
    labels = [row["id"] for row in variants]
    e1_scores = [row.get("e1_score", 0.0) for row in variants]

    colors = ["#1A237E" if e > 0.1 else ("#E65100" if e == 0 else "#616161") for e in e1_scores]
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.bar(x, distances, color=colors, edgecolor="white", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("Composite Bio-Distance (CBD)", fontsize=11)
    ax.set_xlabel("Genomic Variant / Isoform", fontsize=11)
    ax.set_title(
        "Differential Biology Scanner — WT vs Variants (including Locus E1 State)",
        fontsize=12,
        fontweight="bold",
    )
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#1A237E", label="Compartment A"),
        Patch(facecolor="#E65100", label="Boundary Loss (E1 set to 0)"),
        Patch(facecolor="#616161", label="Compartment B")
    ]
    ax.legend(handles=legend_elements, fontsize=9, loc="upper left")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved → {out_png}")


# ── Main ───────────────────────────────────────────────────────────────────────

def main(mode: str = "single"):
    load_dotenv()
    api_key = os.getenv("NVIDIA_ESM_API_KEY") or os.getenv("NVIDIA_API_KEY")
    if not api_key:
        raise EnvironmentError("NVIDIA_ESM_API_KEY (or NVIDIA_API_KEY) not set in .env")

    if mode == "differential":
        if os.path.exists(OUT_MULTI_NPZ) and os.path.exists(OUT_MULTI_JSON):
            print(f"Differential embeddings cache exists → {OUT_MULTI_NPZ}")
            with open(OUT_MULTI_JSON) as f:
                results = json.load(f)
            plot_differential_embeddings(results)
            print("Stage 2 complete.")
            return

        if not os.path.exists(IN_MULTI_FASTA):
            raise FileNotFoundError(f"Multi-FASTA not found: {IN_MULTI_FASTA}. Run Stage 1 differential mode first.")

        records = read_multi_fasta()
        embeddings_by_id = {}
        pooled = {}
        for record in records:
            seq_id = record["header"].split("|", 1)[0]
            sequence = record["sequence"]
            print(f"Calling ESM2-650M for {seq_id} ({len(sequence)} aa) ...")
            embedding = call_esm2_embeddings(sequence, api_key)
            embeddings_by_id[seq_id] = embedding
            pooled[seq_id] = mean_pool_embedding(embedding)

        # ── Cell-State Awareness Integration ──────────────────────────────────────
        context = {}
        if os.path.exists(IN_VARIANT_CONTEXT):
            with open(IN_VARIANT_CONTEXT) as f:
                context = json.load(f)
        
        target_gene = context.get("target_gene", {})
        wt_e1 = target_gene.get("cell_state", {}).get("e1_score", 0.0)
        
        wt = pooled["wt_canonical"]
        results = []
        for record in records:
            seq_id = record["header"].split("|", 1)[0]
            
            # ESM-2 Sequence Distance
            seq_dist = cosine_distance(wt, pooled[seq_id])
            
            # Cell-State Distance (E1 Score mismatch)
            # Logic: If it's the 'boundary_collapse_proxy', it loses its compartment e1 signal.
            if seq_id == "boundary_collapse_proxy":
                v_e1 = 0.0
            else:
                v_e1 = wt_e1
            
            e1_mismatch = abs(wt_e1 - v_e1)
            # Composite Distance (Heuristic: SeqDist + 0.1 * E1 Mismatch)
            # This allows structural genomic 'loss' to amplify the variant score.
            composite_dist = seq_dist + (0.1 * e1_mismatch)
            
            results.append(
                {
                    "id": seq_id,
                    "length": len(record["sequence"]),
                    "embedding_dim": int(pooled[seq_id].shape[0]),
                    "cosine_distance_to_wt": seq_dist,
                    "e1_score": v_e1,
                    "composite_biological_distance": composite_dist
                }
            )

        np.savez(OUT_MULTI_NPZ, **embeddings_by_id)
        with open(OUT_MULTI_JSON, "w") as f:
            json.dump(results, f, indent=2)
        print(f"Differential embeddings saved → {OUT_MULTI_NPZ}")
        print(f"Embedding comparison saved → {OUT_MULTI_JSON}")
        plot_differential_embeddings(results)
        print("Stage 2 complete.")
        return

    if os.path.exists(OUT_NPY):
        print(f"Embeddings cache exists → {OUT_NPY} (skip API call; delete to rerun)")
        try:
            embeddings = np.load(OUT_NPY)
            print(f"  Cached embeddings shape: {embeddings.shape}")
            plot_embeddings(embeddings)
            print("Stage 2 complete.")
            return
        except Exception as exc:
            print(f"  Cache corrupt ({exc}); recomputing...")

    if not os.path.exists(IN_FASTA):
        raise FileNotFoundError(f"FASTA not found: {IN_FASTA}. Run Stage 1 first.")

    sequence = read_fasta()
    print(f"Sequence loaded: {len(sequence)} aa")
    print("Calling ESM2-650M embeddings endpoint...")
    embeddings = call_esm2_embeddings(sequence, api_key)

    os.makedirs(os.path.dirname(OUT_NPY), exist_ok=True)
    np.save(OUT_NPY, embeddings)
    print(f"Embeddings saved → {OUT_NPY}")

    plot_embeddings(embeddings)
    print("Stage 2 complete.")


if __name__ == "__main__":
    main()
