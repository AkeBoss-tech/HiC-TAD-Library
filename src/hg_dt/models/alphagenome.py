import os
from typing import List, Optional

try:
    from alphagenome.models import dna_client
    from alphagenome.data import genome, ontology
except ImportError:
    dna_client = None
    genome = None
    ontology = None


# ---------------------------------------------------------------------------
# Cell-type ontology registry
# Common research cell lines + tissues with their CLO/CL/UBERON terms.
# Add more from https://www.ebi.ac.uk/ols/ontologies/clo
# ---------------------------------------------------------------------------
CELL_TYPES = {
    # Cell lines (CLO)
    "K562":     "CLO:0007050",   # K562 CML cell line
    "GM12878":  "CLO:0002037",   # GM12878 lymphoblastoid
    "HEK293":   "CLO:0001230",   # HEK293 embryonic kidney
    "HeLa":     "CLO:0003684",   # HeLa cervical cancer
    "MCF7":     "CLO:0002756",   # MCF-7 breast cancer
    "Jurkat":   "CLO:0006852",   # Jurkat T-cell leukemia
    "A549":     "CLO:0000149",   # A549 lung carcinoma
    "H1-hESC":  "CLO:0037111",   # H1 human embryonic stem cells
    # Tissues (UBERON)
    "liver":    "UBERON:0002107",
    "brain":    "UBERON:0000955",
    "heart":    "UBERON:0000948",
    "lung":     "UBERON:0002048",
    "kidney":   "UBERON:0002113",
    # Cell types (CL)
    "T-cell":   "CL:0000084",
    "B-cell":   "CL:0000236",
    "neuron":   "CL:0000540",
    "hepatocyte": "CL:0000182",
}

_OUTPUT_MAP = {
    'ATAC':         lambda: dna_client.OutputType.ATAC,
    'CHIP_TF':      lambda: dna_client.OutputType.CHIP_TF,
    'CHIP_HISTONE': lambda: dna_client.OutputType.CHIP_HISTONE,
    'CAGE':         lambda: dna_client.OutputType.CAGE,
    'CONTACT_MAPS': lambda: dna_client.OutputType.CONTACT_MAPS,
    'RNA_SEQ':      lambda: dna_client.OutputType.RNA_SEQ,
    'SPLICE_SITES': lambda: dna_client.OutputType.SPLICE_SITES,
    'SPLICE_SITE_USAGE': lambda: dna_client.OutputType.SPLICE_SITE_USAGE,
    'DNASE':        lambda: dna_client.OutputType.DNASE,
}


def _resolve_outputs(requested_outputs):
    if requested_outputs is None:
        return [
            dna_client.OutputType.ATAC,
            dna_client.OutputType.CHIP_TF,
            dna_client.OutputType.CHIP_HISTONE,
            dna_client.OutputType.CAGE,
            dna_client.OutputType.CONTACT_MAPS,
        ]
    return [_OUTPUT_MAP[r]() for r in requested_outputs if r in _OUTPUT_MAP]


def _resolve_organism(organism: str):
    return dna_client.Organism.HOMO_SAPIENS if organism == "HUMAN" else dna_client.Organism.MUS_MUSCULUS


def _resolve_ontology(cell_type: Optional[str]):
    """Convert a cell type name or CURIE string to an OntologyTerm list."""
    if cell_type is None:
        return None
    if ontology is None:
        return None
    # Allow passing a raw CURIE like "CLO:0007050"
    curie = CELL_TYPES.get(cell_type, cell_type)
    try:
        return [ontology.from_curie(curie)]
    except Exception:
        return None


class AlphaGenomeConnector:
    """
    Wrapper for the AlphaGenome DeepMind API.

    Supports cell-type-specific predictions via ontology terms (CLO, CL, UBERON).
    Use cell_type= with a name from CELL_TYPES or a raw ontology CURIE string.
    """

    def __init__(self, api_key: Optional[str] = None):
        if dna_client is None:
            raise ImportError("AlphaGenome SDK is not installed. Run 'make install-alphagenome'.")

        self.api_key = api_key or os.environ.get("ALPHA_GENOME_API_KEY")
        if not self.api_key:
            raise ValueError("ALPHA_GENOME_API_KEY is not set.")

        self.client = dna_client.create(self.api_key)

    def predict_sequence(
        self,
        sequence: str,
        organism: str = "HUMAN",
        requested_outputs: Optional[List[str]] = None,
        cell_type: Optional[str] = None,
    ):
        """
        Predict functional genomic tracks for a 1 Mb DNA sequence.

        Args:
            sequence: DNA sequence string (~1 Mb)
            organism: "HUMAN" or "MOUSE"
            requested_outputs: subset of ['ATAC', 'CHIP_TF', 'CHIP_HISTONE', 'CAGE',
                               'CONTACT_MAPS', 'RNA_SEQ', 'SPLICE_SITES',
                               'SPLICE_SITE_USAGE', 'DNASE']
            cell_type: name from CELL_TYPES dict (e.g. 'K562') or raw CURIE
                       (e.g. 'CLO:0007050'). None = aggregate across all cell types.
        """
        return self.client.predict_sequence(
            sequence=sequence,
            organism=_resolve_organism(organism),
            requested_outputs=_resolve_outputs(requested_outputs),
            ontology_terms=_resolve_ontology(cell_type),
        )

    def predict_variant(
        self,
        chrom: str,
        start: int,
        end: int,
        variant_pos: int,
        ref_base: str,
        alt_base: str,
        organism: str = "HUMAN",
        requested_outputs: Optional[List[str]] = None,
        cell_type: Optional[str] = None,
    ):
        """
        Predict variant effects for a genomic interval.

        Returns a VariantOutput with .reference and .alternate Output objects.

        Args:
            chrom: chromosome (e.g. 'chr1')
            start/end: 1 Mb window around the variant
            variant_pos: position of the variant (0-based)
            ref_base / alt_base: reference and alternate allele strings
            cell_type: name from CELL_TYPES or raw CURIE (None = all cell types)
        """
        interval = genome.Interval(chromosome=chrom, start=start, end=end)
        variant  = genome.Variant(
            chromosome=chrom,
            position=variant_pos,
            reference_bases=ref_base,
            alternate_bases=alt_base,
        )
        return self.client.predict_variant(
            interval=interval,
            variant=variant,
            organism=_resolve_organism(organism),
            requested_outputs=_resolve_outputs(requested_outputs),
            ontology_terms=_resolve_ontology(cell_type),
        )

    def available_cell_types(self) -> dict:
        """Return the built-in cell type registry."""
        return dict(CELL_TYPES)
