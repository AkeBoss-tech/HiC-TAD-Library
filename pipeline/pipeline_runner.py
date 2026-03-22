#!/usr/bin/env python3
"""
Drug Discovery Pipeline — CLI orchestrator.

Usage:
    python pipeline/pipeline_runner.py            # run all stages
    python pipeline/pipeline_runner.py --stage 2  # single stage
    python pipeline/pipeline_runner.py --force    # ignore caches

Each stage failure is isolated — later stages are still attempted.
"""

import argparse
import importlib
import os
import sys
import time
from datetime import datetime

from dotenv import load_dotenv

# ── Stage registry ─────────────────────────────────────────────────────────────
STAGES = [
    {
        "id": 1,
        "name": "Genomics Context + Sequence Fetch",
        "module": "pipeline.stage1_genomics",
        "outputs": [
            "data/processed/pipeline_target_sequence.fasta",
            "media/pipeline_stage1_expression.png",
        ],
        "api_calls": "0 (UniProt — no key needed)",
        "estimated_time": "~5 s",
    },
    {
        "id": 2,
        "name": "ESM-2 Protein Embeddings",
        "module": "pipeline.stage2_protein",
        "outputs": [
            "data/processed/pipeline_esm2_embeddings.npy",
            "media/pipeline_stage2_embeddings.png",
        ],
        "api_calls": "1 (NVIDIA NIM ESM2-30M)",
        "estimated_time": "~20 s",
    },
    {
        "id": 3,
        "name": "ESMFold Structure Prediction",
        "module": "pipeline.stage3_folding",
        "outputs": [
            "data/processed/pipeline_structure.pdb",
            "media/pipeline_stage3_structure.png",
        ],
        "api_calls": "1 (NVIDIA NIM ESMFold)",
        "estimated_time": "~60 s",
    },
    {
        "id": 4,
        "name": "Molecule Generation + Docking",
        "module": "pipeline.stage4_molecules",
        "outputs": [
            "data/processed/pipeline_molecules.json",
            "media/pipeline_stage4_molecules.png",
        ],
        "api_calls": "22 (1 MolMIM + 1 ESM2 + 20 DiffDock)",
        "estimated_time": "~10 min",
    },
]


# ── Helpers ────────────────────────────────────────────────────────────────────

def check_env(stages_to_run: list[dict]) -> None:
    """Check required API keys for the requested stages."""
    load_dotenv()
    ids = [s["id"] for s in stages_to_run]
    if ids == [1]:
        return  # Stage 1 only uses UniProt (no key)

    # Per-model keys with fallback to NVIDIA_API_KEY
    fallback = os.getenv("NVIDIA_API_KEY")
    checks = {
        2: ("NVIDIA_ESM_API_KEY",      os.getenv("NVIDIA_ESM_API_KEY")      or fallback),
        3: ("NVIDIA_ESM_FOLD_API_KEY", os.getenv("NVIDIA_ESM_FOLD_API_KEY") or fallback),
        4: ("NVIDIA_MOLMIM_API_KEY",   os.getenv("NVIDIA_MOLMIM_API_KEY")   or fallback),
    }
    missing = []
    for stage_id, (var, val) in checks.items():
        if stage_id in ids:
            status = f"✓ (len={len(val)})" if val else "❌ MISSING"
            print(f"  {var}: {status}")
            if not val:
                missing.append(var)

    if missing:
        raise EnvironmentError(
            f"Missing API keys: {missing}. "
            "Set per-model keys or a single NVIDIA_API_KEY in .env."
        )


def clear_caches(stages_to_run: list[dict]) -> None:
    """Delete output files so stages rerun from scratch."""
    for stage in stages_to_run:
        for path in stage["outputs"]:
            if os.path.exists(path):
                os.remove(path)
                print(f"  Removed cache: {path}")


def run_stage(stage: dict) -> tuple[str, float]:
    """Import and run stage main(). Returns (status, elapsed_seconds)."""
    t0 = time.time()
    try:
        mod = importlib.import_module(stage["module"])
        importlib.reload(mod)
        mod.main()
        status = "✅ OK"
    except Exception as exc:
        status = f"❌ FAILED: {exc}"
    elapsed = time.time() - t0
    return status, elapsed


def print_summary(results: list[dict]) -> None:
    """Print ASCII summary table."""
    print("\n" + "=" * 72)
    print(f"{'Pipeline Summary':^72}")
    print("=" * 72)
    fmt = "{:<6} {:<34} {:<14} {}"
    print(fmt.format("Stage", "Name", "Time", "Status"))
    print("-" * 72)
    for r in results:
        elapsed_str = f"{r['elapsed']:.1f}s"
        name_trunc  = r["name"][:33]
        print(fmt.format(f"[{r['id']}]", name_trunc, elapsed_str, r["status"]))
    print("=" * 72)

    failed = [r for r in results if "FAILED" in r["status"]]
    if failed:
        print(f"\n{len(failed)} stage(s) failed. Check output above for details.")
    else:
        print("\nAll stages completed successfully.")

    total = sum(r["elapsed"] for r in results)
    print(f"Total elapsed: {total:.1f}s\n")


# ── CLI ────────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="Drug Discovery Pipeline — UNC5B death domain",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python pipeline/pipeline_runner.py            # all stages
  python pipeline/pipeline_runner.py --stage 2  # only Stage 2
  python pipeline/pipeline_runner.py --force    # ignore caches
        """,
    )
    parser.add_argument(
        "--stage", type=int, choices=[1, 2, 3, 4], default=None,
        help="Run a single stage (1–4). Omit to run all stages.",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Delete existing output caches before running.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Add repo root to sys.path so pipeline.* imports work
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)

    stages_to_run = STAGES if args.stage is None else [s for s in STAGES if s["id"] == args.stage]

    print(f"\nDrug Discovery Pipeline — {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Stages: {[s['id'] for s in stages_to_run]}\n")

    check_env(stages_to_run)

    if args.force:
        print("--force: clearing caches ...")
        clear_caches(stages_to_run)

    results = []
    for stage in stages_to_run:
        print(f"\n{'─' * 60}")
        print(f"[Stage {stage['id']}] {stage['name']}")
        print(f"  API calls       : {stage['api_calls']}")
        print(f"  Estimated time  : {stage['estimated_time']}")
        print(f"{'─' * 60}")
        status, elapsed = run_stage(stage)
        results.append({"id": stage["id"], "name": stage["name"],
                         "status": status, "elapsed": elapsed})

    print_summary(results)


if __name__ == "__main__":
    main()
