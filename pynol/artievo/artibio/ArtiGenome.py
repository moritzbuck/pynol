from pynol.artievo.artibio.ArtiGene import ArtiGene
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

class ArtiGenome(object):
    mutate = lambda seq, r : "".join([s if uniform(0,1) > r else choice(list(set(BASES).difference(s))) for s in seq])

    def __get__(self, i):
        return self.genome[i]

    def __iter__(self): return iter(self.genome)
    def __getitem__(self, key): return self.genome[key]

    @classmethod
    def simplest(cls, genome_length):
        return ArtiGenome(1, genome_length, 1)

    def simi(self, other):
        assert len(self.genome) == len(other.genome), "only works for genomes of same length right now"
        assert all([len(g2.seq) == len(g2.seq) for g1, g2 in zip(self.genome, other.genome)]), "only works for genomes of same length right now"
        a = self.full_seq()
        b = other.full_seq()
        return sum([x == y for x,y in zip(a,b)]) / len(a)


    def __init__(self, nb_genes, length_genes, functional_cutoffs):
        self.nb_genes = nb_genes
        self.length_genes = length_genes if type(length_genes) == list else [length_genes] * self.nb_genes
        self.functional_cutoffs = functional_cutoffs if type(functional_cutoffs) == list else [functional_cutoffs] * self.nb_genes

        self.genome = [ArtiGene.random(m,f) for i,m,f in zip(range(self.nb_genes),self.length_genes, self.functional_cutoffs)]

    def copy(self):
        out_genome = ArtiGenome(self.nb_genes,self.length_genes, self.functional_cutoffs)
        out_genome.genome = [ArtiGene.copy(g) for g in self.genome]
        return out_genome

    def full_seq(self):
        return "".join([s.seq for s in self.genome])

    def to_seqrec(self):
        return SeqRecord(Seq(self.full_seq()), description = "", id = "ArtiGenome_{id:X}".format(id = id(self)))
