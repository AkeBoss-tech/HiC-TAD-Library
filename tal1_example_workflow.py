#!/usr/bin/env python3
"""
AlphaGenome TAL1 locus example workflow.

Implements the workflow from:
https://www.alphagenomedocs.com/colabs/example_analysis_workflow.html

- Load TAL1 oncogenic variants and visualize their positions
- Predict effect of one variant (Jurkat) on RNA-seq, DNase, ChIP-histone
- Compare predicted TAL1 expression: cancer-associated vs shuffled background variants
"""

import io
import itertools
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from dotenv import load_dotenv

# Ensure alphagenome is on path when run from project root (e.g. pip install -e externals/alphagenome)
from alphagenome.data import gene_annotation, genome
from alphagenome.data import transcript as transcript_utils
from alphagenome.models import dna_client
from alphagenome.models import variant_scorers
from alphagenome.visualization import plot_components

load_dotenv()
API_KEY = os.getenv("ALPHA_GENOME_API_KEY")
if not API_KEY:
    try:
        from alphagenome import colab_utils
        API_KEY = colab_utils.get_api_key()
    except Exception:
        raise ValueError("ALPHA_GENOME_API_KEY not found. Set it in .env or environment.") from None

# Output directory
MEDIA_DIR = "media"
os.makedirs(MEDIA_DIR, exist_ok=True)


# ---------- Data loading ----------

def load_gtf_and_extractor():
    """Load GENCODE GTF and build MANE-select transcript extractor."""
    gtf_url = (
        "https://storage.googleapis.com/alphagenome/reference/gencode/"
        "hg38/gencode.v46.annotation.gtf.gz.feather"
    )
    print("Loading gene annotations (GENCODE hg38)...")
    gtf = pd.read_feather(gtf_url)
    gtf_transcript = gene_annotation.filter_protein_coding(gtf)
    gtf_transcript = gene_annotation.filter_to_mane_select_transcript(gtf_transcript)
    extractor = transcript_utils.TranscriptExtractor(gtf_transcript)
    return extractor


# ---------- Variant data (TAL1 oncogenic) ----------

def oncogenic_tal1_variants() -> pd.DataFrame:
    """Returns a dataframe of oncogenic T-ALL variants that affect TAL1."""
    variant_data = """
ID	CHROM	POS	REF	ALT	output	Study ID	Study Variant ID
Jurkat	chr1	47239296	C	CCGTTTCCTAACC	1	Mansour_2014
MOLT-3	chr1	47239296	C	ACC	1	Mansour_2014
Patient_1	chr1	47239296	C	AACG	1	Mansour_2014
Patient_2	chr1	47239291	CTAACC	TTTACCGTCTGTTAACGGC	1	Mansour_2014
Patient_3-5	chr1	47239296	C	ACG	1	Mansour_2014
Patient_6	chr1	47239296	C	ACC	1	Mansour_2014
Patient_7	chr1	47239295	AC	TCAAACTGGTAACC	1	Mansour_2014
Patient_8	chr1	47239296	C	AACC	1	Mansour_2014
new 3' enhancer 1	chr1	47212072	T	TGGGTAAACCGTCTGTTCAGCG	1	Smith_2023	UPNT802
new 3' enhancer 2	chr1	47212074	G	GAACGTT	1	Smith_2023	UPNT613
intergenic SNV 1	chr1	47230639	C	T	1	Liu_2020	SJALL043861_D1
intergenic SNV 2	chr1	47230639	C	T	1	Liu_2020	SJALL018373_D1
SJALL040467_D1	chr1	47239296	C	AACC	1	Liu_2020	SJALL040467_D1
PATBGC	chr1	47239296	C	AACC	1	Liu_2017	PATBGC
PATBTX	chr1	47239296	C	ACGGATATAACC	1	Liu_2017	PATBTX
PARJAY	chr1	47239296	C	ACGGAATTTCTAACC	1	Liu_2017	PARJAY
PARSJG	chr1	47239296	C	AACC	1	Liu_2017	PARSJG
PASYAJ	chr1	47239296	C	AACC	1	Liu_2017	PASYAJ
PATRAB	chr1	47239293	TTA	CTAACGG	1	Liu_2017	PATRAB
PAUBXP	chr1	47239296	C	ACC	1	Liu_2017	PAUBXP
PATENL	chr1	47239296	C	AACC	1	Liu_2017	PATENL
PARNXJ	chr1	47239296	C	ACG	1	Liu_2017	PARNXJ
PASXSI	chr1	47239296	C	AACC	1	Liu_2017	PASXSI
PASNEH	chr1	47239296	C	ACC	1	Liu_2017	PASNEH
PAUAFN	chr1	47239296	C	AACC	1	Liu_2017	PAUAFN
PARASZ	chr1	47239296	C	ACC	1	Liu_2017	PARASZ
PARWNW	chr1	47239296	C	ACC	1	Liu_2017	PARWNW
PASFKA	chr1	47239293	TTA	ACCGTTAATCAA	1	Liu_2017	PASFKA
PATEIT	chr1	47239296	C	AC	1	Liu_2017	PATEIT
PASMHF	chr1	47239296	C	AC	1	Liu_2017	PASMHF
PARJNX	chr1	47239296	C	AC	1	Liu_2017	PARJNX
PASYWF	chr1	47239296	C	AC	1	Liu_2017	PASYWF
"""
    return pd.read_table(io.StringIO(variant_data.strip()), sep="\t")


def vcf_row_to_variant(row: pd.Series) -> genome.Variant:
    """Parse a row of a VCF-style dataframe into genome.Variant."""
    return genome.Variant(
        chromosome=str(row["CHROM"]),
        position=int(row["POS"]),
        reference_bases=row["REF"],
        alternate_bases=row["ALT"],
        name=row["ID"],
    )


def generate_background_variants(variant: genome.Variant, max_number: int = 100) -> pd.DataFrame:
    """Generate background (shuffled) variants of same length at same position."""
    nucleotides = np.array(list("ACGT"), dtype="<U1")

    def generate_unique_strings(n, max_n, random_seed=42):
        rng = np.random.default_rng(random_seed)
        if 4**n < max_n:
            raise ValueError("Cannot generate that many unique strings for the given length.")
        generated = set()
        while len(generated) < max_n:
            indices = rng.integers(0, 4, size=n)
            new_string = "".join(nucleotides[indices])
            if new_string != variant.alternate_bases:
                generated.add(new_string)
        return list(generated)

    if 4 ** len(variant.alternate_bases) < max_number:
        permutations = ["".join(p) for p in itertools.product(nucleotides, repeat=len(variant.alternate_bases))]
    else:
        permutations = generate_unique_strings(len(variant.alternate_bases), max_number)

    return pd.DataFrame({
        "ID": ["mut_" + str(variant.position) + "_" + x for x in permutations],
        "CHROM": variant.chromosome,
        "POS": variant.position,
        "REF": variant.reference_bases,
        "ALT": permutations,
        "output": 0.0,
        "original_variant": variant.name,
    })


def inference_df(qtl_df: pd.DataFrame, input_sequence_length: int) -> pd.DataFrame:
    """Build dataframe of variants and intervals for inference."""
    rows = []
    for _, row in qtl_df.iterrows():
        variant = vcf_row_to_variant(row)
        interval = genome.Interval(
            chromosome=row["CHROM"], start=row["POS"], end=row["POS"]
        ).resize(input_sequence_length)
        rows.append({
            "interval": interval,
            "variant": variant,
            "output": row["output"],
            "variant_id": row["ID"],
            "POS": row["POS"],
            "REF": row["REF"],
            "ALT": row["ALT"],
            "CHROM": row["CHROM"],
        })
    return pd.DataFrame(rows)


def oncogenic_and_background_variants(
    input_sequence_length: int, number_of_background_variants: int = 20
) -> pd.DataFrame:
    """Oncogenic + shuffled background variants for evaluation."""
    oncogenic = oncogenic_tal1_variants()
    variants = [
        genome.Variant(
            chromosome=str(r.CHROM),
            position=int(r.POS),
            reference_bases=r.REF,
            alternate_bases=r.ALT,
            name=r.ID,
        )
        for r in oncogenic.itertuples()
    ]
    background = pd.concat([
        generate_background_variants(v, number_of_background_variants) for v in variants
    ])
    all_df = pd.concat([oncogenic, background])
    return inference_df(all_df, input_sequence_length=input_sequence_length)


def coarse_grained_mute_groups(eval_df: pd.DataFrame):
    """Group variants by position and ALT length for plotting."""
    grp = []
    for row in eval_df.itertuples():
        if row.POS >= 47239290:  # MUTE site
            grp.append("MUTE_other" if row.ALT_len > 4 else "MUTE_" + str(row.ALT_len))
        else:
            grp.append(str(row.POS) + "_" + str(row.ALT_len))
    return pd.Categorical(grp, categories=sorted(set(grp)), ordered=True)


# ---------- Main workflow ----------

def main():
    print("AlphaGenome TAL1 example workflow")
    print("=" * 50)

    dna_model = dna_client.create(API_KEY)
    transcript_extractor = load_gtf_and_extractor()

    # TAL1 interval (~30 kb around TAL1)
    tal1_interval = genome.Interval(
        chromosome="chr1", start=47209255, end=47242023, strand="-"
    )

    # 1) Visualize variant positions
    print("\n1. Plotting variant positions near TAL1...")
    unique_positions = np.sort(oncogenic_tal1_variants()["POS"].unique())
    labels = [
        "47212072, 47212074", "", "47230639", "47239291 - 47239296", "", "", ""
    ][: len(unique_positions)]
    _ = plot_components.plot(
        [
            plot_components.TranscriptAnnotation(transcript_extractor.extract(tal1_interval)),
        ],
        annotations=[
            plot_components.VariantAnnotation(
                [
                    genome.Variant(
                        chromosome="chr1", position=int(x), reference_bases="N", alternate_bases="N"
                    )
                    for x in unique_positions
                ],
                labels=labels[: len(unique_positions)],
                use_default_labels=False,
            )
        ],
        interval=tal1_interval,
        title="Positions of variants near TAL1",
    )
    out1 = os.path.join(MEDIA_DIR, "tal1_variant_positions.png")
    plt.savefig(out1, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"   Saved: {out1}")

    # 2) Predict effect of one variant (Jurkat) on RNA-seq, DNase, ChIP-histone
    print("\n2. Predicting effect of Jurkat variant (RNA-seq, DNase, ChIP-histone)...")
    variant = vcf_row_to_variant(oncogenic_tal1_variants().iloc[0])
    ontology_terms = ["CL:0001059"]  # CD34+ common myeloid progenitor

    output = dna_model.predict_variant(
        interval=tal1_interval.resize(2**20),
        variant=variant,
        requested_outputs={
            dna_client.OutputType.RNA_SEQ,
            dna_client.OutputType.CHIP_HISTONE,
            dna_client.OutputType.DNASE,
        },
        ontology_terms=ontology_terms,
    )

    transcripts = transcript_extractor.extract(tal1_interval)
    _ = plot_components.plot(
        [
            plot_components.TranscriptAnnotation(transcripts),
            plot_components.Tracks(
                tdata=output.alternate.rna_seq.filter_to_nonpositive_strand()
                - output.reference.rna_seq.filter_to_nonpositive_strand(),
                ylabel_template="{biosample_name} ({strand})\n{name}",
                filled=True,
            ),
            plot_components.Tracks(
                tdata=output.alternate.dnase.filter_to_nonpositive_strand()
                - output.reference.dnase.filter_to_nonpositive_strand(),
                ylabel_template="{biosample_name} ({strand})\n{name}",
                filled=True,
            ),
            plot_components.Tracks(
                tdata=output.alternate.chip_histone.filter_to_nonpositive_strand()
                - output.reference.chip_histone.filter_to_nonpositive_strand(),
                ylabel_template="{biosample_name} ({strand})\n{name}",
                filled=True,
            ),
        ],
        annotations=[plot_components.VariantAnnotation([variant])],
        interval=tal1_interval,
        title=(
            "Effect of variant on predicted RNA Expression, DNAse, and ChIP-Histone "
            f"in CD34+ HSC.\n{variant=}"
        ),
    )
    out2 = os.path.join(MEDIA_DIR, "tal1_jurkat_variant_effect.png")
    plt.savefig(out2, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"   Saved: {out2}")

    # 3) Score variants: oncogenic vs background (small set for speed)
    print("\n3. Scoring variants (oncogenic vs background)...")
    n_background = 2  # keep small for a quick run
    eval_df = oncogenic_and_background_variants(
        input_sequence_length=2**20, number_of_background_variants=n_background
    )
    eval_df["ALT_len"] = eval_df["ALT"].str.len()
    eval_df["variant_group"] = eval_df["POS"].astype(str) + "_" + eval_df["ALT_len"].astype(str)
    eval_df["output"] = eval_df["output"].fillna(0) != 0
    eval_df["coarse_grained_variant_group"] = coarse_grained_mute_groups(eval_df)

    scores = dna_model.score_variants(
        intervals=eval_df["interval"].to_list(),
        variants=eval_df["variant"].to_list(),
        variant_scorers=[variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]],
        max_workers=2,
    )

    gene_index = scores[0][0].obs.query('gene_name == "TAL1"').index[0]
    cell_type_index = scores[0][0].var.query('ontology_curie == "CL:0001059"').index[0]

    def get_tal1_score(score_data):
        return score_data[gene_index, cell_type_index].X[0, 0]

    eval_df["tal1_diff_in_cd34"] = [get_tal1_score(x[0]) for x in scores]
    print("   TAL1 expression change (CD34+): computed for all variants.")

    # 4) Comparison plot (plotnine if available, else simple matplotlib)
    try:
        import plotnine as gg
    except ImportError:
        gg = None

    plot_df = eval_df.loc[eval_df.REF != eval_df.ALT].copy()
    plot_df["variant"] = plot_df["variant"].astype(str)
    plot_df = plot_df[["variant", "output", "tal1_diff_in_cd34", "coarse_grained_variant_group"]].drop_duplicates()

    if gg is not None and len(plot_df.coarse_grained_variant_group.unique()) > 0:
        print("\n4. Building comparison plots (cancer vs background)...")
        facet_title_by_group = {
            "47212072_22": "chr1:47212072\n21 bp ins.",
            "47212074_7": "chr1:47212072\n7 bp",
            "47230639_1": "chr1:47230639\nSNV",
            "MUTE_2": "chr1:47239296\n1 bp ins.",
            "MUTE_3": "chr1:47239296\n2 bp ins.",
            "MUTE_4": "chr1:47239296\n3 bp ins.",
            "MUTE_other": "chr1:47239296\n7-18 bp ins.",
        }
        for group in plot_df.coarse_grained_variant_group.unique():
            sub = plot_df[plot_df.coarse_grained_variant_group == group]
            if len(sub) == 0:
                continue
            sub = pd.concat([
                sub.assign(plot_group="density"),
                sub.assign(plot_group="rain"),
            ])
            sub = sub[~((sub.plot_group == "density") & sub.output)]
            col_width = max(np.ptp(sub.tal1_diff_in_cd34) / 200, 1e-6)
            sub["col_width"] = sub["output"].map({True: 1.5 * col_width, False: 1.25 * col_width})
            title = facet_title_by_group.get(str(group), str(group))
            p = (
                gg.ggplot(sub)
                + gg.aes(x="tal1_diff_in_cd34")
                + gg.geom_col(
                    gg.aes(y=1, width="col_width", fill="output", x="tal1_diff_in_cd34", alpha="output"),
                    data=sub[sub["plot_group"] == "rain"],
                )
                + gg.geom_density(
                    gg.aes(x="tal1_diff_in_cd34", fill="output"),
                    data=sub[sub["plot_group"] == "density"],
                    color="white",
                )
                + gg.facet_wrap("~output + plot_group", nrow=1, scales="free_x")
                + gg.scale_alpha_manual({True: 1, False: 0.3})
                + gg.scale_fill_manual({True: "#FAA41A", False: "gray"})
                + gg.labs(title=title)
                + gg.theme_minimal()
                + gg.geom_vline(xintercept=0, linetype="dotted")
                + gg.theme(
                    figure_size=(1.2, 3),
                    legend_position="none",
                    axis_text_x=gg.element_blank(),
                    panel_grid_major_x=gg.element_blank(),
                    panel_grid_minor_x=gg.element_blank(),
                    strip_text=gg.element_blank(),
                    axis_title_y=gg.element_blank(),
                    axis_title_x=gg.element_blank(),
                    plot_title=gg.element_text(size=9),
                )
                + gg.scale_y_reverse()
                + gg.coord_flip()
            )
            safe_name = str(group).replace(" ", "_")
            out3 = os.path.join(MEDIA_DIR, f"tal1_comparison_{safe_name}.png")
            p.save(out3, dpi=150, width=4, height=3)
            print(f"   Saved: {out3}")
    else:
        if gg is None:
            print("\n4. Skipping plotnine comparison (install plotnine for comparison plots).")
        # Save a simple summary table
        summary_path = os.path.join(MEDIA_DIR, "tal1_scores_summary.csv")
        plot_df.to_csv(summary_path, index=False)
        print(f"   Saved score summary: {summary_path}")

    print("\nDone. Check the 'media/' directory for outputs.")


if __name__ == "__main__":
    main()
    sys.exit(0)
