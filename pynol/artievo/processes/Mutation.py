from random import random
from random import choices
from pynol.artievo.artibio import BASES
from pynol.artievo.artibio.ArtiGene import ArtiGene
from pynol.artievo.artibio.ArtiGenome import ArtiGenome

default_probs = { (a,b) : 1.0/(len(BASES)-1) for a in BASES for b in BASES}


class Mutation:
    def __init__(self, rate):
        self.rate = rate

    def mutate(gene):
        pass

class PointMutation(Mutation):

    def __init__(self, rate, substitution_matrix = default_probs):
        Mutation.__init__(self, rate)
        self.substitution_matrix = substitution_matrix

    def _substitute(self, base):
        possibles = list(set(BASES).difference(base))
        new_base = choices(possibles, [self.substitution_matrix[(base, p)] for p in possibles] )
        return new_base[0]

    def mutate_gene(self, gene):
        out_gene = ArtiGene.copy(gene)
        out_gene.seq = "".join([s if random() > self.rate else self._substitute(s) for s in out_gene.seq])
        return out_gene

    def mutate_genome(self, genome):
        out_genome = ArtiGenome.copy(genome)
        out_genome.genome = [self.mutate_gene(g) for g in out_genome]
        return out_genome
