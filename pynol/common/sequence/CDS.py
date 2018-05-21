# File: Genomic.py

from pynol.common.sequence.Nucleotide import Nucleotide

class CDS( Nucleotide ):

    def get_prot(self):
        self.sequence.translate()
