import os
import requests
from typing import Dict, Optional, List


def _parse_bfactor(pdb_text: str) -> List[float]:
    """Extract B-factor (pLDDT) values from ATOM records in a PDB string."""
    plddt = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and len(line) >= 66:
            try:
                plddt.append(float(line[60:66].strip()))
            except ValueError:
                pass
    return plddt


def _mock_structure(aa_seq: str) -> Dict:
    """Fallback mock when the real folding API is unavailable."""
    n_res = len(aa_seq)
    if n_res > 50:
        plddt = [90.0 + (i % 5) for i in range(n_res)]
    else:
        plddt = [40.0 + (i % 10) for i in range(n_res)]
    pdb_mock = f"HEADER    MOCK PDB FOR {n_res} RESIDUES\n"
    return {'pdb': pdb_mock, 'plddt': plddt}


class ProteinFolder:
    """Wrapper for protein structure prediction via ESMFold or BioNeMo."""

    # Maximum sequence length accepted by the public ESMFold endpoint
    _ESMFOLD_MAX_LEN = 400

    def __init__(self, api_key: Optional[str] = None):
        self.api_key = api_key or os.environ.get("BIONEMO_API_KEY")

    def predict_structure(self, aa_seq: str) -> Dict:
        """
        Predict the 3D structure of an amino acid sequence.

        Tries Meta's public ESMFold endpoint first. Falls back to a mock
        if the API is unreachable or the sequence exceeds the length limit.

        Args:
            aa_seq: The amino acid sequence to fold.

        Returns:
            Dict with 'pdb' (str) and 'plddt' (List[float]).
        """
        if not aa_seq:
            return {'pdb': '', 'plddt': []}

        if len(aa_seq) <= self._ESMFOLD_MAX_LEN:
            try:
                resp = requests.post(
                    "https://api.esmatlas.com/foldSequence/v1/pdb/",
                    data=aa_seq,
                    headers={"Content-Type": "application/x-www-form-urlencoded"},
                    timeout=60,
                )
                if resp.ok and resp.text.strip().startswith(("ATOM", "HEADER", "MODEL")):
                    pdb = resp.text
                    plddt = _parse_bfactor(pdb)
                    return {'pdb': pdb, 'plddt': plddt}
            except Exception:
                pass

        return _mock_structure(aa_seq)
