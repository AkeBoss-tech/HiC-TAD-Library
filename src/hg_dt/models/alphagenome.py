import os
from typing import List, Optional, Iterable
import numpy as np

try:
    from alphagenome.models import dna_client
    from alphagenome.data import genome
except ImportError:
    # Handle optional installation
    dna_client = None
    genome = None


class AlphaGenomeConnector:
    """Wrapper for the AlphaGenome DeepMind API."""

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
        organism="HUMAN",
        requested_outputs: Optional[List[str]] = None
    ):
        """
        Predict all relevant tracks for a 1 Mb sequence.

        Args:
            sequence: 1 Mb DNA sequence
            organism: "HUMAN" or "MOUSE"
            requested_outputs: list of strings (e.g. ['ATAC', 'CHIP_TF', 'CHIP_HISTONE', 'CAGE', 'CONTACT_MAPS'])
        """
        if requested_outputs is None:
            requested_outputs = [
                dna_client.OutputType.ATAC,
                dna_client.OutputType.CHIP_TF,
                dna_client.OutputType.CHIP_HISTONE,
                dna_client.OutputType.CAGE,
                dna_client.OutputType.CONTACT_MAPS
            ]
        else:
            # Map string to enum
            output_map = {
                'ATAC': dna_client.OutputType.ATAC,
                'CHIP_TF': dna_client.OutputType.CHIP_TF,
                'CHIP_HISTONE': dna_client.OutputType.CHIP_HISTONE,
                'CAGE': dna_client.OutputType.CAGE,
                'CONTACT_MAPS': dna_client.OutputType.CONTACT_MAPS,
                'RNA_SEQ': dna_client.OutputType.RNA_SEQ
            }
            requested_outputs = [output_map[req] for req in requested_outputs if req in output_map]

        org = dna_client.Organism.HOMO_SAPIENS if organism == "HUMAN" else dna_client.Organism.MUS_MUSCULUS

        output = self.client.predict_sequence(
            sequence=sequence,
            organism=org,
            requested_outputs=requested_outputs,
            ontology_terms=None # Predict all if None
        )
        return output

    def predict_variant(
        self,
        chrom: str,
        start: int,
        end: int,
        variant_pos: int,
        ref_base: str,
        alt_base: str,
        organism="HUMAN",
        requested_outputs: Optional[List[str]] = None
    ):
        """Predict variant effects for a 1 Mb sequence."""
        if requested_outputs is None:
            requested_outputs = [
                dna_client.OutputType.ATAC,
                dna_client.OutputType.CHIP_TF,
                dna_client.OutputType.CHIP_HISTONE,
                dna_client.OutputType.CAGE,
                dna_client.OutputType.CONTACT_MAPS
            ]
        else:
            # Map string to enum
            output_map = {
                'ATAC': dna_client.OutputType.ATAC,
                'CHIP_TF': dna_client.OutputType.CHIP_TF,
                'CHIP_HISTONE': dna_client.OutputType.CHIP_HISTONE,
                'CAGE': dna_client.OutputType.CAGE,
                'CONTACT_MAPS': dna_client.OutputType.CONTACT_MAPS,
                'RNA_SEQ': dna_client.OutputType.RNA_SEQ
            }
            requested_outputs = [output_map[req] for req in requested_outputs if req in output_map]

        interval = genome.Interval(chromosome=chrom, start=start, end=end)
        variant = genome.Variant(
            chromosome=chrom,
            position=variant_pos,
            reference_bases=ref_base,
            alternate_bases=alt_base
        )

        org = dna_client.Organism.HOMO_SAPIENS if organism == "HUMAN" else dna_client.Organism.MUS_MUSCULUS

        outputs = self.client.predict_variant(
            interval=interval,
            variant=variant,
            organism=org,
            requested_outputs=requested_outputs,
            ontology_terms=None
        )
        return outputs
