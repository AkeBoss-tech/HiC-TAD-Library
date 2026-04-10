"""
UCSC Genome Browser REST API client.

Fetches hg38 sequences and gene annotations without requiring local FASTA/GTF files.
All data comes from the public UCSC REST API — no API key required.

API docs: https://api.genome.ucsc.edu/
"""

import re
import requests
import pandas as pd
from functools import lru_cache
from typing import Optional, Tuple

_UCSC_BASE = "https://api.genome.ucsc.edu"
_TIMEOUT   = 30  # seconds

# AlphaGenome sequence input must be exactly this many bp (2^20).
ALPHA_WINDOW_SIZE = 1_048_576


def build_alpha_window(
    center_bp: int,
    size: int = ALPHA_WINDOW_SIZE,
    chrom_len: Optional[int] = None,
) -> Tuple[int, int]:
    """
    Return [start, end) genomic interval (0-based, half-open) centered on
    ``center_bp``, with length ``size`` (default 1 Mb for AlphaGenome).

    If ``chrom_len`` is given, the window is clamped so ``end <= chrom_len``.
    Otherwise behavior matches the previous app (only ``start`` is floored at 0).
    """
    start = max(0, center_bp - size // 2)
    end = start + size
    if chrom_len is not None:
        end = min(end, chrom_len)
        start = max(0, end - size)
    return start, end


def normalize_insert_sequence(raw: str, expected_len: int) -> str:
    """
    Uppercase ACGT only; strip whitespace. If empty, return empty string
    (caller may substitute a placeholder).
    """
    s = re.sub(r"\s+", "", (raw or "").upper())
    if not s:
        return ""
    if not re.fullmatch(r"[ACGT]+", s):
        raise ValueError("Insert sequence must contain only A/C/G/T.")
    if len(s) != expected_len:
        raise ValueError(
            f"Insert sequence length {len(s)} bp must match edit span {expected_len} bp."
        )
    return s


def compute_mut_shift(mod_type: str, edit_start: int, edit_end: int, insert_len: int) -> int:
    """Net bp change in mutant vs reference for coordinate mapping (see ``translate.transcript``)."""
    t = (mod_type or "").lower()
    if t == "deletion":
        return -(edit_end - edit_start)
    if t == "insertion":
        return insert_len
    return 0


def fetch_ref_mut_sequences(
    chrom: str,
    edit_start: int,
    edit_end: int,
    mod_type: str,
    *,
    insert_seq: str = "",
    snp_alt: Optional[str] = None,
    window_size: int = ALPHA_WINDOW_SIZE,
    fasta_path: Optional[str] = None,
) -> Tuple[str, str, int, int]:
    """
    Fetch hg38 windows and build ref + mut sequences, both exactly ``window_size`` bp.

    Parameters
    ----------
    insert_seq
        For *insertion*, the inserted bases (length must equal ``edit_end - edit_start``),
        or empty to use a run of ``A`` of that length.
    snp_alt
        For *snp*, alternate base ``A``/``C``/``G``/``T``. If None, picks the first
        base different from the reference at the SNP position (legacy behavior).
    fasta_path
        If set, read sequence from this hg38 FASTA (:func:`fetch_hg38_sequence_local`)
        instead of the UCSC REST API.
    """
    edit_size = edit_end - edit_start
    if edit_size <= 0 and mod_type.lower() != "none":
        raise ValueError("edit_end must be greater than edit_start.")

    center = (edit_start + edit_end) // 2
    win_start, win_end = build_alpha_window(center, size=window_size)

    def pull(c: str, s: int, e: int) -> str:
        if fasta_path:
            return fetch_hg38_sequence_local(c, s, e, fasta_path)
        return fetch_hg38_sequence(c, s, e)

    ref_seq = pull(chrom, win_start, win_end)
    t = mod_type.lower()

    if t == "deletion":
        mut_raw = pull(chrom, win_start, win_end + edit_size)
        rel_start = edit_start - win_start
        rel_end = edit_end - win_start
        mut_seq = mut_raw[:rel_start] + mut_raw[rel_end:]

    elif t == "insertion":
        ins = normalize_insert_sequence(insert_seq, edit_size) if insert_seq else ""
        if not ins:
            ins = "A" * edit_size
        fetch_end = max(win_start + 1000, win_end - edit_size)
        mut_raw = pull(chrom, win_start, fetch_end)
        rel_start = edit_start - win_start
        mut_seq = mut_raw[:rel_start] + ins + mut_raw[rel_start:]

    elif t == "snp":
        rel_start = edit_start - win_start
        orig = ref_seq[rel_start]
        if snp_alt:
            alt = snp_alt.upper()
            if alt not in "ACGT":
                raise ValueError("SNP alternate must be A, C, G, or T.")
            if alt == orig.upper():
                raise ValueError("SNP alternate must differ from the reference base.")
        else:
            alt = next(b for b in ("A", "T", "C", "G") if b != orig.upper())
        mut_seq = ref_seq[:rel_start] + alt + ref_seq[rel_start + 1 :]

    else:
        mut_seq = ref_seq

    if len(ref_seq) != window_size:
        raise RuntimeError(f"ref_seq length {len(ref_seq)} != {window_size}")
    if len(mut_seq) != window_size:
        raise RuntimeError(f"mut_seq length {len(mut_seq)} != {window_size}")
    return ref_seq, mut_seq, win_start, win_end


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
    # AlphaGenome max window is 1,048,576 (2^20); allow extra for deletion compensation
    if end - start > 2_000_000:
        raise ValueError(
            f"Requested region too large: {end - start:,} bp. Max 2,000,000 bp."
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


@lru_cache(maxsize=8)
def _open_local_fasta(path: str):
    from pyfaidx import Fasta

    return Fasta(path)


def _resolve_fasta_chrom(fa, chrom: str) -> Optional[str]:
    """Match UCSC-style ``chr1`` names to pyfaidx keys."""
    for c in (chrom, chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"):
        if c in fa:
            return c
    return None


def fetch_hg38_sequence_local(chrom: str, start: int, end: int, fasta_path: str) -> str:
    """
    Read hg38 sequence from a local FASTA (e.g. UCSC ``hg38.fa``) using pyfaidx.
    Coordinates are 0-based half-open, matching :func:`fetch_hg38_sequence`.
    """
    if end - start > 2_000_000:
        raise ValueError(
            f"Requested region too large: {end - start:,} bp. Max 2,000,000 bp."
        )
    fa = _open_local_fasta(fasta_path)
    ck = _resolve_fasta_chrom(fa, chrom)
    if ck is None:
        raise RuntimeError(f"Chromosome {chrom!r} not found in {fasta_path}")
    chrom_len = len(fa[ck])
    if end > chrom_len:
        raise ValueError(
            f"End {end} past chromosome {ck} length {chrom_len}"
        )
    return str(fa[ck][start:end]).upper()


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


@lru_cache(maxsize=64)
def fetch_gene_models(chrom: str, start: int, end: int) -> list:
    """
    Fetch detailed gene models (with exon coordinates) from UCSC ncbiRefSeq.

    Returns a list of dicts, each with:
      name, name2, strand, txStart, txEnd, cdsStart, cdsEnd,
      exon_starts (list[int]), exon_ends (list[int])
    """
    url = (
        f"{_UCSC_BASE}/getData/track"
        f"?genome=hg38&track=ncbiRefSeq&chrom={chrom}&start={start}&end={end}"
    )
    try:
        resp = requests.get(url, timeout=_TIMEOUT)
        if not resp.ok:
            return []
        records = resp.json().get("ncbiRefSeq", [])
    except Exception:
        return []

    models = []
    seen = set()
    for r in records:
        gene  = r.get("name2") or r.get("name", "")
        tx_id = r.get("name", "")
        key   = (tx_id, r.get("txStart"), r.get("txEnd"))
        if key in seen:
            continue
        seen.add(key)

        # Parse comma-separated exon lists
        def _parse_list(s):
            return [int(x) for x in str(s).split(",") if x.strip()]

        models.append({
            "name":        tx_id,
            "name2":       gene,
            "strand":      r.get("strand", "+"),
            "txStart":     int(r.get("txStart", 0)),
            "txEnd":       int(r.get("txEnd",   0)),
            "cdsStart":    int(r.get("cdsStart", r.get("txStart", 0))),
            "cdsEnd":      int(r.get("cdsEnd",   r.get("txEnd",   0))),
            "exon_starts": _parse_list(r.get("exonStarts", "")),
            "exon_ends":   _parse_list(r.get("exonEnds",   "")),
        })
    return models


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
