#!/usr/bin/env python3
"""Summarize pairtools dedup stats for a LinkPrep preprocessing run."""

from __future__ import annotations

import argparse
import gzip
import json
from dataclasses import dataclass, asdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate LinkPrep QC summaries.")
    parser.add_argument("--stats", required=True, help="Path to pairtools dedup stats file")
    parser.add_argument("--pairs", required=False, help="Path to mapped.pairs.gz")
    parser.add_argument("--mcool", required=False, help="Path to generated .mcool")
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument("--markdown-out", required=True, help="Output markdown path")
    parser.add_argument("--json-out", required=True, help="Output JSON path")
    return parser.parse_args()


def parse_stats(stats_path: Path) -> dict[str, int]:
    metrics: dict[str, int] = {}
    for line in stats_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        key, value = parts[0], parts[1]
        try:
            metrics[key] = int(value)
        except ValueError:
            continue
    return metrics


def count_pairs(pairs_path: Path | None) -> int | None:
    if pairs_path is None or not pairs_path.exists():
        return None
    count = 0
    with gzip.open(pairs_path, "rt") as handle:
        for line in handle:
            if not line.startswith("#"):
                count += 1
    return count


@dataclass
class QCSummary:
    sample: str
    stats_path: str
    pairs_path: str | None
    mcool_path: str | None
    total_pairs: int | None
    no_dup_pairs: int | None
    duplicate_pairs: int | None
    cis_1kb_pairs: int | None
    pairs_file_records: int | None
    duplicate_rate_pct: float | None
    cis_1kb_of_no_dup_pct: float | None
    passes_cis_1kb_threshold: bool | None
    passes_depth_threshold: bool | None


def build_summary(args: argparse.Namespace) -> QCSummary:
    stats_path = Path(args.stats)
    pairs_path = Path(args.pairs) if args.pairs else None
    mcool_path = Path(args.mcool) if args.mcool else None
    metrics = parse_stats(stats_path)

    total_pairs = metrics.get("total")
    no_dup_pairs = metrics.get("no_dup")
    duplicate_pairs = metrics.get("duplicate")
    cis_1kb_pairs = metrics.get("cis_1kb+")
    if cis_1kb_pairs is None:
        cis_1kb_pairs = metrics.get("cis_1kb")

    duplicate_rate_pct = None
    if total_pairs and duplicate_pairs is not None:
        duplicate_rate_pct = duplicate_pairs * 100.0 / total_pairs

    cis_1kb_of_no_dup_pct = None
    if no_dup_pairs and cis_1kb_pairs is not None:
        cis_1kb_of_no_dup_pct = cis_1kb_pairs * 100.0 / no_dup_pairs

    passes_cis_1kb_threshold = None if cis_1kb_of_no_dup_pct is None else cis_1kb_of_no_dup_pct > 40.0
    passes_depth_threshold = None if no_dup_pairs is None else no_dup_pairs > 125_000_000

    return QCSummary(
        sample=args.sample,
        stats_path=str(stats_path),
        pairs_path=str(pairs_path) if pairs_path else None,
        mcool_path=str(mcool_path) if mcool_path else None,
        total_pairs=total_pairs,
        no_dup_pairs=no_dup_pairs,
        duplicate_pairs=duplicate_pairs,
        cis_1kb_pairs=cis_1kb_pairs,
        pairs_file_records=count_pairs(pairs_path),
        duplicate_rate_pct=duplicate_rate_pct,
        cis_1kb_of_no_dup_pct=cis_1kb_of_no_dup_pct,
        passes_cis_1kb_threshold=passes_cis_1kb_threshold,
        passes_depth_threshold=passes_depth_threshold,
    )


def format_pct(value: float | None) -> str:
    return "n/a" if value is None else f"{value:.2f}%"


def status_label(value: bool | None) -> str:
    if value is None:
        return "n/a"
    return "PASS" if value else "CHECK"


def write_markdown(summary: QCSummary, out_path: Path) -> None:
    out_path.write_text(
        "\n".join(
            [
                f"# LinkPrep QC Summary: {summary.sample}",
                "",
                "## Outputs",
                "",
                f"- Stats: `{summary.stats_path}`",
                f"- Pairs: `{summary.pairs_path or 'n/a'}`",
                f"- mcool: `{summary.mcool_path or 'n/a'}`",
                "",
                "## Metrics",
                "",
                f"- Total pairs: `{summary.total_pairs}`",
                f"- No-dup pairs: `{summary.no_dup_pairs}`",
                f"- Duplicate pairs: `{summary.duplicate_pairs}`",
                f"- cis >= 1 kb pairs: `{summary.cis_1kb_pairs}`",
                f"- Records in `.pairs.gz`: `{summary.pairs_file_records}`",
                f"- Duplicate rate: `{format_pct(summary.duplicate_rate_pct)}`",
                f"- cis >= 1 kb as % of no-dup: `{format_pct(summary.cis_1kb_of_no_dup_pct)}`",
                "",
                "## Threshold Checks",
                "",
                f"- Dovetail ligation threshold (`cis >= 1 kb > 40% of no-dup`): `{status_label(summary.passes_cis_1kb_threshold)}`",
                f"- Dovetail depth threshold (`no-dup > 125M`): `{status_label(summary.passes_depth_threshold)}`",
                "",
                "## Interpretation",
                "",
                "- `PASS` means the observed metric exceeds the rule-of-thumb threshold.",
                "- `CHECK` means the run completed, but the sample should be reviewed before trusting downstream TAD/loop sensitivity.",
                "- `n/a` means the required metric was not present in the stats file.",
                "",
            ]
        )
    )


def write_json(summary: QCSummary, out_path: Path) -> None:
    out_path.write_text(json.dumps(asdict(summary), indent=2))


def main() -> None:
    args = parse_args()
    summary = build_summary(args)
    write_markdown(summary, Path(args.markdown_out))
    write_json(summary, Path(args.json_out))


if __name__ == "__main__":
    main()
