import os
from pyfaidx import Fasta
import pyranges as pr
import pandas as pd

class ReferenceContextBuilder:
    def __init__(self, fasta_path, gtf_path, bed_path):
        """
        Initialize the ReferenceContextBuilder.

        Args:
            fasta_path (str): Path to the reference FASTA file or directory.
            gtf_path (str): Path to the GENCODE GTF file.
            bed_path (str): Path to the ENCODE cCREs BED file.
        """
        self.fasta_path = fasta_path
        import glob

        # Determine if fasta_path is a directory
        if os.path.isdir(fasta_path):
            fa_files = glob.glob(os.path.join(fasta_path, "*.fa")) + glob.glob(os.path.join(fasta_path, "*.fa.gz"))
            if not fa_files:
                raise FileNotFoundError(f"No FASTA files found in directory {fasta_path}")
            # Fasta class from pyfaidx accepts a list of fasta files in recent versions,
            # or we can construct a unified dictionary of Fasta objects.
            # Pyfaidx's Fasta class natively accepts a single filename string.
            # If multiple are passed, we have to create our own multiplexer.

            # Simple Dictionary Multiplexer:
            class MultiFasta:
                def __init__(self, files):
                    self.fastas = [Fasta(f) for f in files]
                    self.chroms = {}
                    for f in self.fastas:
                        for chrom in f.keys():
                            self.chroms[chrom] = f[chrom]
                def __getitem__(self, key):
                    return self.chroms[key]

            self.fasta = MultiFasta(fa_files)
        else:
            self.fasta = Fasta(fasta_path)

        # Load annotations using pyranges
        self.gtf = pr.read_gtf(gtf_path)

        # Bed files may lack headers or standard column names depending on format
        try:
            self.bed = pr.read_bed(bed_path)
        except Exception:
            # Fallback for some custom BEDs
            df = pd.read_csv(bed_path, sep='\t', header=None, comment='#')
            df.columns = ['Chromosome', 'Start', 'End'] + [f'Col_{i}' for i in range(3, len(df.columns))]
            self.bed = pr.PyRanges(df)

    def _get_window(self, chrom, start, end, window_size=1000000):
        """
        Calculate a window centered on the edit.
        """
        center = (start + end) // 2
        half_window = window_size // 2

        # Ensure we don't go below 0 (0-based indexing)
        win_start = max(0, center - half_window)
        win_end = win_start + window_size

        return win_start, win_end

    def _get_mutant_seq(self, ref_seq, edit_start, edit_end, edit_type, new_sequence, win_start):
        """
        Apply the edit to the reference sequence string.
        """
        # Convert absolute coordinates to 0-based relative coordinates within the window
        rel_start = edit_start - win_start
        rel_end = edit_end - win_start

        if edit_type == 'deletion':
            mut_seq = ref_seq[:rel_start] + ref_seq[rel_end:]
        elif edit_type == 'insertion':
            mut_seq = ref_seq[:rel_start] + new_sequence + ref_seq[rel_start:]
        elif edit_type == 'SNP' or edit_type == 'replacement':
            mut_seq = ref_seq[:rel_start] + new_sequence + ref_seq[rel_end:]
        else:
            raise ValueError(f"Unknown edit_type: {edit_type}")

        return mut_seq

    def get_context(self, chrom, start, end, edit_type, new_sequence=None, window_size=1000000):
        """
        Get 1 Mb context around an edit.

        Args:
            chrom (str): Chromosome (e.g., 'chr17')
            start (int): Start coordinate of the edit (0-based)
            end (int): End coordinate of the edit (0-based, exclusive)
            edit_type (str): 'SNP', 'deletion', 'insertion', 'replacement'
            new_sequence (str): The sequence to insert or replace with (optional for deletions)
            window_size (int): Size of the surrounding window (default 1 Mb)

        Returns:
            dict containing:
                - ref_seq: The reference FASTA string for the window
                - mut_seq: The mutant FASTA string for the window
                - annotations: Dictionary of overlapping 'genes' and 'cCREs' DataFrames
        """
        win_start, win_end = self._get_window(chrom, start, end, window_size)

        # 1. Get Reference Sequence
        # pyfaidx uses 0-based, half-open intervals [start, end)
        ref_seq_obj = self.fasta[chrom][win_start:win_end]
        ref_seq = str(ref_seq_obj)

        # 2. Generate Mutant Sequence
        mut_seq = self._get_mutant_seq(ref_seq, start, end, edit_type, new_sequence, win_start)

        # 3. Get Annotations
        # Create a single-interval PyRanges object for the window
        window_pr = pr.from_dict({"Chromosome": [chrom], "Start": [win_start], "End": [win_end]})

        # Find overlapping GTF features
        overlapping_gtf = self.gtf.intersect(window_pr)
        # Find overlapping BED features
        overlapping_bed = self.bed.intersect(window_pr)

        # Extract features as pandas DataFrames (can be empty)
        gtf_df = overlapping_gtf.df if not overlapping_gtf.empty else pd.DataFrame()
        bed_df = overlapping_bed.df if not overlapping_bed.empty else pd.DataFrame()

        # Try to subset standard relevant features if present
        genes_df = pd.DataFrame()
        if not gtf_df.empty and 'Feature' in gtf_df.columns:
            genes_df = gtf_df[gtf_df['Feature'].isin(['exon', 'gene', 'transcript'])]
        elif not gtf_df.empty:
            genes_df = gtf_df

        return {
            'ref_seq': ref_seq,
            'mut_seq': mut_seq,
            'annotations': {
                'genes': genes_df,
                'cCREs': bed_df
            },
            'window': (win_start, win_end)
        }
