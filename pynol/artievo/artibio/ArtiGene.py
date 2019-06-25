from random import choices
from pynol.artievo.artibio import BASES

class ArtiGene:
    def __init__(self, seq, functional_cutoff):
        self.seq = seq
        self.functional_cutoff = functional_cutoff

    @classmethod
    def copy(cls, gene):
        outp = cls(gene.seq, gene.functional_cutoff)
        return outp

    @classmethod
    def random(cls, length, functional_cutoff):
        seq = "".join(choices(BASES, k=length))
        return cls(seq, functional_cutoff)
