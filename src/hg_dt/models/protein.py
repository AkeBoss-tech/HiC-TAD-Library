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
    """Fallback mock when all folding APIs are unavailable."""
    n_res = len(aa_seq)
    plddt = [90.0 + (i % 5) for i in range(n_res)] if n_res > 50 else [40.0 + (i % 10) for i in range(n_res)]
    return {'pdb': f"HEADER    MOCK PDB FOR {n_res} RESIDUES\n", 'plddt': plddt, 'source': 'mock'}


class ProteinFolder:
    """
    Protein structure prediction via NVIDIA ESMFold NIM (primary) or Meta public
    ESMFold endpoint (secondary), with mock fallback.

    Priority:
      1. NVIDIA ESMFold NIM  — fastest, reliable, uses NVIDIA_ESM_FOLD_API_KEY
      2. Meta public ESMFold — free, up to 400 AA, no key needed
      3. Mock                — synthetic pLDDT based on length
    """

    _NVIDIA_ENDPOINT = "https://health.api.nvidia.com/v1/biology/nvidia/esmfold"
    _META_ENDPOINT   = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    _META_MAX_LEN    = 400  # Meta public ESMFold input limit
    ESMFOLD_MAX_INPUT_AA = _META_MAX_LEN  # alias for dashboards / warnings

    def __init__(self, api_key: Optional[str] = None):
        self.api_key = (
            api_key
            or os.environ.get("NVIDIA_ESM_FOLD_API_KEY")
            or os.environ.get("NVIDIA_ESM_API_KEY")
            or os.environ.get("NVIDIA_API_KEY")
        )

    def predict_structure(self, aa_seq: str) -> Dict:
        """
        Predict the 3D structure of an amino acid sequence.

        Returns dict with:
          'pdb'    : PDB file string
          'plddt'  : list of per-residue confidence scores (0-100)
          'source' : which backend was used
        """
        if not aa_seq:
            return {'pdb': '', 'plddt': [], 'source': 'none'}

        # 1. Try NVIDIA ESMFold NIM
        if self.api_key:
            try:
                resp = requests.post(
                    self._NVIDIA_ENDPOINT,
                    json={"sequence": aa_seq},
                    headers={"Authorization": f"Bearer {self.api_key}",
                             "Content-Type": "application/json"},
                    timeout=120,
                )
                if resp.ok:
                    data = resp.json()
                    pdb = data["pdbs"][0] if data.get("pdbs") else ""
                    if pdb:
                        return {'pdb': pdb, 'plddt': _parse_bfactor(pdb), 'source': 'nvidia_esmfold'}
            except Exception:
                pass

        # 2. Try Meta public ESMFold
        if len(aa_seq) <= self._META_MAX_LEN:
            try:
                resp = requests.post(
                    self._META_ENDPOINT,
                    data=aa_seq,
                    headers={"Content-Type": "application/x-www-form-urlencoded"},
                    timeout=60,
                )
                if resp.ok and resp.text.strip().startswith(("ATOM", "HEADER", "MODEL", "PARENT")):
                    pdb = resp.text
                    return {'pdb': pdb, 'plddt': _parse_bfactor(pdb), 'source': 'meta_esmfold'}
            except Exception:
                pass

        # 3. Mock fallback
        return _mock_structure(aa_seq)
