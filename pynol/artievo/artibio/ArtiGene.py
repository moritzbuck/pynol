from random import choices
from pynol.artievo.artibio import BASES

class ArtiGene:
    def __init__(self, seq):
        self.seq = "".join(choices(BASES, k=length_gene))

    @classmethod
    def from_gene(cls, gene):
        return cls(gene.seq)

    @classmethod
    def random(cls, length):
        seq = "".join(choices(BASES, k=length))
        return cls(seq)
