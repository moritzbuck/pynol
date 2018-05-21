# File: Nucleotide.py

from pynol.common.sequence.Sequence import Sequence
from Bio.Alphabet import DNAAlphabet

class Nucleotide(Sequence):

# Attributes: Instance
# Operations

    @Sequence.sequence.setter
    def sequence(self, value):
        assert value.alphabet.__class__ == DNAAlphabet, "This is not a nucleotide sequence"
        Sequence.sequence.fset(self,value)
