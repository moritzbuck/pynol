import os
from os.path import join as pjoin
from Bio.Alphabet import DNAAlphabet
from Bio.Seq import Seq
from pynol.common.sequence.Genomic import Genomic
from pynol.common.genome.Genome import Genome
from pynol.common.genome.sources.FromRefSeq import FromRefSeq
from pynol.common.genome.sources.FromFile import FromFile
from pynol.common.taxonomy.Taxon import Taxon
from pynol.common.taxonomy.Taxonomy import Taxonomy
from pynol.database.PynolDB import PynolDB
from pynol.common.sequence.RNA import RNA
from pynol.common.sequence.CDS import CDS
from pynol.common.sequence.Feature import Feature
import pynol
from tqdm import tqdm
from pynol.tools.gene_prediction.Prokka import Prokka
from pynol.common.COGs.COGing import COGing
from pynol.common.COGs.COG import COG
from bson.objectid import ObjectId
import pandas


from mongo_thingy import connect, Thingy
connect("mongodb://localhost/pateci")

taxonomy = Taxonomy()
taxonomy.from_gtdb_taxonomy_file("/home/moritzbuck/data/bac_taxonomy_r80.tsv")

prok = Prokka.find_one()

all_gtdb_patesc = COGing.find_one()

search_string = "^{tax_string}.*".format(tax_string = taxonomy['p__Patescibacteria'].get_tax_string(full = True))
all_pates = list(Genome.find({'taxonomy.gtdb' : {'$regex' : search_string}}))

all_pate_taxa = [t for t in taxonomy.taxa_dict.values() if t.is_child(taxonomy['p__Patescibacteria'])]
all_pate = {tt : list(Genome.find({'taxonomy.gtdb' : {'$regex' : tt.get_tax_string(full=True) + '.*'}})) for tt in all_pate_taxa}


files = [n for n in  os.listdir("test_data/heatmap_and_cores/") if ".cogs" in n]

cores = {}
for f in files :
     with open("test_data/heatmap_and_cores/" + f) as handle:
         cores[f.split(".")[0]] = [l.strip() for l in handle]

all_cores = set(sum(list(cores.values()),[]))

core_matrix = { c : {g : 1 for g in COG.find_one(ObjectId(c)).genomes} for c in all_cores}
core_matrix = pandas.DataFrame.from_dict(core_matrix, orient="index").fillna(0)
core_matrix.to_csv("test_data/core_matrix.csv")

cog_metadata = {k : {kk : int(kk in cc)  for kk in core_matrix.index} for k, cc in cores.items()}
cog_metadata = pandas.DataFrame.from_dict(cog_metadata)
cog_metadata.to_csv("test_data/cog_metadata.csv")

genome_taxonomy = pandas.DataFrame.from_dict({ g.id : {Taxon.rev_leves[i] : taxonomy[g.taxonomy['gtdb']].get_parent(i, "full") for i in range(3,6) if taxonomy[g.taxonomy['gtdb']].get_parent(i, "full") in cores.keys() }   for g in all_pates}, orient = "index").fillna("other")
genome_taxonomy.to_csv("test_data/genome_taxonomy.csv")

genome_md = pandas.DataFrame.from_dict({g.id : { t : (taxonomy[g.taxonomy['gtdb']].is_child(taxonomy[t]))+1-1   for t in cores.keys() } for g in all_pates }, orient = "index").fillna(0)
genome_md['genome_length'] = [Genome.find_one(ObjectId(g)).corrected_length for g in genome_md.index]
genome_md.to_csv("test_data/genome_md.csv")
with open("test_data/genome_taxonomy.csv","w") as handle :
        handle.writelines([",taxon\n"])
        handle.writelines(["{id},{tax}\n".format(id = g.id,  tax = g.taxonomy['gtdb']) for g in all_pates])


#for taxa, genomes  in all_pate.items():
#    if len(genomes) > 50:
#        print(taxa)
#        matrix = all_gtdb_patesc.coocurrences(cutoff=len(genomes) * 0.1, genome_set = set([g.id for g in genomes]))
#        with open("test_data/" + taxa.get_tax_string(full =True).split(";")[-1] + ".clstrs", "w") as handle:
#            handle.writelines(["COG_1,COG_2,overlap\n"])
#            handle.writelines(["{g1},{g2},{dist}\n".format(g1=kk[0], g2=kk[1], dist=vv) for kk,vv in  matrix.items() ])
