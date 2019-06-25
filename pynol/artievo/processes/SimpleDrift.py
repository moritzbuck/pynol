from pynol.artievo.artibio.ArtiGenome import ArtiGenome
from pynol.artievo.processes.Mutation import PointMutation
from tqdm import tqdm
from random import random
from Bio import SeqIO

class SimpleDrift:

    def __init__(self, genome_length = 1000, mutation_rate = 0.001, death_rate = 0.01, birth_rate = 0.012, start_size = 100):
        self.seed_genome = ArtiGenome.simplest(genome_length)
        self.seed_genome.origin = (None, None)
        self.mutation = PointMutation(mutation_rate)
        self.population = [self.mutation.mutate_genome(self.seed_genome) for i in range(start_size)]
        self.alls = [p.copy() for p in self.population]
        self.gen = 0
        self.death_rate = death_rate
        self.birth_rate = birth_rate

    def run(self, nb_gens = 400):
        for i in range(nb_gens):
            self.one_gen()
            self.gen += 1

    def one_gen(self):
        babies = []
        kill_bill = []
        for g in self.population:
            if random() < self.birth_rate:
                baby = self.mutation.mutate_genome(g)
                baby.origin = (g, g.simi(baby))
                babies += [baby]
            if random() <  self.death_rate:
                kill_bill += [g]
        for g in kill_bill:
            self.population.remove(g)
        self.population += babies
        self.alls += [b.copy() for b in babies]

        avg_simi = sum([self.seed_genome.simi(g) for g in self.population])/len(self.population)
        print("Generation {i}, nb_genomes : {count}, {simi}".format(i = self.gen, count =len(self.population),simi = avg_simi))
