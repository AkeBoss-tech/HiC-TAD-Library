"""
UCSC Genome Browser REST API client.

Fetches hg38 sequences and gene annotations without requiring local FASTA/GTF files.
All data comes from the public UCSC REST API — no API key required.

API docs: https://api.genome.ucsc.edu/
"""

import requests
import pandas as pd
from functools import lru_cache

_UCSC_BASE = "https://api.genome.ucsc.edu"
_TIMEOUT   = 30  # seconds


@lru_cache(maxsize=128)
def fetch_hg38_sequence(chrom: str, start: int, end: int) -> str:
    """
    Fetch the hg38 reference DNA sequence for the given interval.

    Parameters
    ----------
    chrom : str   e.g. "chr1"
    start : int   0-based start (UCSC convention)
    end   : int   0-based end (exclusive)

    Returns
    -------
    str  — uppercase DNA sequence (A/T/C/G/N)

    Raises
    ------
    RuntimeError if the API call fails.
    """
    # AlphaGenome requires power-of-2 lengths; max supported is 1,048,576 (2^20)
    if end - start > 1_100_000:
        raise ValueError(
            f"Requested region too large: {end - start:,} bp. "
            "Use 1,048,576 bp (AlphaGenome max window)."
        )

    url = (
        f"{_UCSC_BASE}/getData/sequence"
        f"?genome=hg38&chrom={chrom}&start={start}&end={end}"
    )
    resp = requests.get(url, timeout=_TIMEOUT)
    if not resp.ok:
        raise RuntimeError(
            f"UCSC sequence fetch failed: HTTP {resp.status_code} — {resp.text[:200]}"
        )
    data = resp.json()
    seq = data.get("dna", "")
    if not seq:
        raise RuntimeError(f"UCSC returned empty sequence for {chrom}:{start}-{end}")
    return seq.upper()


@lru_cache(maxsize=64)
def fetch_gene_annotations(chrom: str, start: int, end: int) -> pd.DataFrame:
    """
    Fetch NCBI RefSeq gene annotations in the given interval from UCSC.

    Returns a DataFrame with columns: Name, Type, Start, End.
    The Start/End are genomic coordinates (0-based).
    Returns an empty DataFrame if no genes are found.
    """
    url = (
        f"{_UCSC_BASE}/getData/track"
        f"?genome=hg38&track=ncbiRefSeq&chrom={chrom}&start={start}&end={end}"
    )
    try:
        resp = requests.get(url, timeout=_TIMEOUT)
        if not resp.ok:
            return pd.DataFrame(columns=["Name", "Type", "Start", "End"])
        data = resp.json()
    except Exception:
        return pd.DataFrame(columns=["Name", "Type", "Start", "End"])

    # The response has a key matching the track name (e.g. "ncbiRefSeq")
    records = data.get("ncbiRefSeq", [])
    if not records:
        return pd.DataFrame(columns=["Name", "Type", "Start", "End"])

    rows = []
    seen = set()
    for r in records:
        name = r.get("name2") or r.get("name", "unknown")
        tx_start = r.get("txStart", 0)
        tx_end   = r.get("txEnd", 0)
        key = (name, tx_start, tx_end)
        if key in seen:
            continue
        seen.add(key)
        rows.append({"Name": name, "Type": "Gene", "Start": tx_start, "End": tx_end})

    if not rows:
        return pd.DataFrame(columns=["Name", "Type", "Start", "End"])
    return pd.DataFrame(rows)


@lru_cache(maxsize=64)
def fetch_regulatory_elements(chrom: str, start: int, end: int) -> pd.DataFrame:
    """
    Fetch ENCODE cCREs (candidate cis-Regulatory Elements) from the UCSC
    encodeCcreCombined track for hg38.

    Returns DataFrame with columns: Name, Type, Start, End.
    Type is one of: 'Enhancer', 'Promoter', 'CTCF', 'Insulator', 'cCRE'.
    """
    url = (
        f"{_UCSC_BASE}/getData/track"
        f"?genome=hg38&track=encodeCcreCombined&chrom={chrom}&start={start}&end={end}"
    )
    try:
        resp = requests.get(url, timeout=_TIMEOUT)
        if not resp.ok:
            return pd.DataFrame(columns=["Name", "Type", "Start", "End"])
        data = resp.json()
    except Exception:
        return pd.DataFrame(columns=["Name", "Type", "Start", "End"])

    records = data.get("encodeCcreCombined", [])
    if not records:
        return pd.DataFrame(columns=["Name", "Type", "Start", "End"])

    # Map ENCODE cCRE descriptions to readable types
    _type_map = {
        "PLS":   "Promoter",
        "pELS":  "Enhancer",
        "dELS":  "Enhancer",
        "CTCF-only": "CTCF",
        "DNase-H3K4me3": "Promoter",
    }

    rows = []
    for i, r in enumerate(records[:50]):   # cap at 50 elements for display
        ccre_type = r.get("ucscLabel") or r.get("description", "cCRE")
        readable = _type_map.get(ccre_type, "cCRE")
        rows.append({
            "Name": f"{readable}_{i+1}",
            "Type": readable,
            "Start": r.get("chromStart", r.get("start", start)),
            "End":   r.get("chromEnd",  r.get("end",   end)),
        })

    if not rows:
        return pd.DataFrame(columns=["Name", "Type", "Start", "End"])
    return pd.DataFrame(rows)


def build_locus_annotations(chrom: str, edit_start: int, edit_end: int,
                             window: int = 500_000) -> pd.DataFrame:
    """
    Build a combined annotation DataFrame for Step 2 preview.
    Fetches real genes + cCREs from UCSC around the edit locus.

    Falls back to an empty DataFrame if both calls fail (e.g. offline).
    """
    pad = window // 2
    view_start = max(0, (edit_start + edit_end) // 2 - pad)
    view_end   = (edit_start + edit_end) // 2 + pad

    genes = fetch_gene_annotations(chrom, view_start, view_end)
    ccres = fetch_regulatory_elements(chrom, view_start, view_end)

    parts = [df for df in [genes, ccres] if not df.empty]
    if not parts:
        return pd.DataFrame(columns=["Name", "Type", "Start", "End"])
    return pd.concat(parts, ignore_index=True)
