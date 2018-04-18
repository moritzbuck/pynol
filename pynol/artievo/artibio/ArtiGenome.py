from pynol.artievo.artibio import BASES

class ArtiGenome(object):
    mutate = lambda seq, r : "".join([s if uniform(0,1) > r else choice(list(set(BASES).difference(s))) for s in seq])
    simi = lambda a,b : sum([x == y for x,y in zip(a,b)])/len(a)

    def __get__(self, i):
        return self.genome[i]
    def __iter__(self): return iter(self.genome)
    def __getitem__(self, key): return self.genome[key]

    def __init__(self, nb_genes, length_genes, model = None):
        self.nb_genes = nb_genes
        self.length_genes = length_genes
        self.model = model
        self.genome = [Gene(self.length_genes, self.mutation_rates, self.functional_cutoffs) for i,m,f in zip(range(self.nb_genes), self.mutation_rates, self.functional_cutoffs)]

    def copy(self):
        out_genome = Genome(self.nb_genes,self.length_genes)
        out_genome.genome = self.genome
