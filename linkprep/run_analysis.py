"""
Master entry point for LinkPrep data analysis.

Steps:
  1. Generate sample data (if needed)
  2. QC metrics from pairs file
  3. Contact analysis: TADs, compartments
  4. SV detection
  5. CNV analysis
  6. Haplotype phasing

Usage:
    # With sample data (generated automatically):
    python run_analysis.py

    # With real data (update config.py first):
    python run_analysis.py --no-generate

    # Run only specific steps:
    python run_analysis.py --steps qc contacts sv cnv phasing
"""

import argparse
import os
import subprocess
import sys
import time


def banner(title: str):
    width = 60
    print("\n" + "=" * width)
    print(f"  {title}")
    print("=" * width)


def run_step(name: str, fn, *args, **kwargs):
    banner(name)
    t0 = time.time()
    try:
        result = fn(*args, **kwargs)
        print(f"  Done in {time.time() - t0:.1f}s")
        return result
    except Exception as exc:
        print(f"  ERROR: {exc}")
        import traceback
        traceback.print_exc()
        return None


def main():
    parser = argparse.ArgumentParser(description="LinkPrep analysis pipeline")
    parser.add_argument("--no-generate", action="store_true",
                        help="Skip sample data generation (use real data)")
    parser.add_argument("--steps", nargs="+",
                        choices=["generate", "qc", "contacts", "sv", "cnv", "phasing"],
                        default=["generate", "qc", "contacts", "sv", "cnv", "phasing"],
                        help="Which analysis steps to run")
    args = parser.parse_args()

    # Ensure we can import from the linkprep package
    here = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.dirname(here)
    sys.path.insert(0, repo)
    sys.path.insert(0, here)

    import config

    # ------------------------------------------------------------------
    # Step 0: Generate sample data
    # ------------------------------------------------------------------
    if "generate" in args.steps and not args.no_generate:
        banner("Generating sample data")
        gen_script = os.path.join(here, "sample_data", "generate_sample_data.py")
        if not os.path.exists(config.COOL_FILE):
            t0 = time.time()
            subprocess.run([sys.executable, gen_script], check=True)
            print(f"  Done in {time.time() - t0:.1f}s")
        else:
            print("  Sample data already exists — skipping generation.")
            print(f"  (Delete {config.DATA_DIR} to regenerate)")

    # ------------------------------------------------------------------
    # Check required files
    # ------------------------------------------------------------------
    missing = []
    for label, path in [
        ("Cool file",    config.COOL_FILE),
        ("Pairs file",   config.PAIRS_FILE),
        ("VCF file",     config.VCF_FILE),
        ("Coverage BED", config.COVERAGE_BED),
    ]:
        if not os.path.exists(path):
            missing.append(f"  {label}: {path}")
    if missing:
        print("\nMissing required files:")
        print("\n".join(missing))
        print("\nRun without --no-generate, or update config.py with real data paths.")
        sys.exit(1)

    # ------------------------------------------------------------------
    # Step 1: QC
    # ------------------------------------------------------------------
    if "qc" in args.steps:
        from pipeline.qc_metrics import read_pairs_sample, compute_qc, print_qc
        from pipeline.qc_metrics import plot_distance_distribution, plot_pair_type_bar

        def qc_step():
            df    = read_pairs_sample(config.PAIRS_FILE)
            stats = compute_qc(df)
            print_qc(stats)
            plot_distance_distribution(
                stats["cis_distances"],
                os.path.join(config.FIGURES_DIR, "qc_distance_distribution.png"),
            )
            plot_pair_type_bar(
                stats,
                os.path.join(config.FIGURES_DIR, "qc_pair_types.png"),
            )

        run_step("QC Metrics", qc_step)

    # ------------------------------------------------------------------
    # Step 2: Contact analysis
    # ------------------------------------------------------------------
    if "contacts" in args.steps:
        import cooler
        from visualization.plot_contacts import run_all as run_contacts

        def contacts_step():
            clr = cooler.Cooler(config.COOL_FILE)
            run_contacts(clr)

        run_step("Contact Analysis (TADs + Compartments)", contacts_step)

    # ------------------------------------------------------------------
    # Step 3: SV detection
    # ------------------------------------------------------------------
    if "sv" in args.steps:
        from visualization.plot_sv import run_all as run_sv

        def sv_step():
            run_sv(config.PAIRS_FILE)

        run_step("Structural Variant Detection", sv_step)

    # ------------------------------------------------------------------
    # Step 4: CNV
    # ------------------------------------------------------------------
    if "cnv" in args.steps:
        from visualization.plot_cnv import run_all as run_cnv
        run_step("Copy-Number Variation", run_cnv)

    # ------------------------------------------------------------------
    # Step 5: Phasing
    # ------------------------------------------------------------------
    if "phasing" in args.steps:
        from visualization.plot_phasing import run_all as run_phasing
        run_step("Haplotype Phasing", run_phasing)

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    banner("Analysis Complete")
    print(f"  All figures saved to: {config.FIGURES_DIR}")
    print()
    figures = sorted(os.listdir(config.FIGURES_DIR))
    for f in figures:
        if f.endswith(".png"):
            print(f"    {f}")
    print()


if __name__ == "__main__":
    main()
