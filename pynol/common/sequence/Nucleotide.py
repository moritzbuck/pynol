# File: Nucleotide.py

from pynol.common.sequence.Sequence import Sequence
from Bio.Alphabet import DNAAlphabet

class Nucleotide(Sequence):

# Attributes: Instance

    def __init__(self):
        super(Sequence, self).__init__()
# Operations

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, value):
        assert value.alphabet.__class__ == DNAAlphabet, "This is not a nucleotide sequence"
        self._sequence = value
