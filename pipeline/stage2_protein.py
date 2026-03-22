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
OUT_NPY  = "data/processed/pipeline_esm2_embeddings.npy"
OUT_PNG  = "media/pipeline_stage2_embeddings.png"


# ── Helpers ────────────────────────────────────────────────────────────────────

def read_fasta(path: str = IN_FASTA) -> str:
    lines = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith(">"):
                lines.append(line)
    return "".join(lines)


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


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    load_dotenv()
    api_key = os.getenv("NVIDIA_ESM_API_KEY") or os.getenv("NVIDIA_API_KEY")
    if not api_key:
        raise EnvironmentError("NVIDIA_ESM_API_KEY (or NVIDIA_API_KEY) not set in .env")

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
