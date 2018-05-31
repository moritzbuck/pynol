import os
from os.path import join as pjoin
from pynol.common.sequence.Genomic import Genomic
from pynol.common.genome.Genome import Genome
from pynol.common.taxonomy.Taxonomy import Taxonomy
from tqdm import tqdm

from mongo_thingy import connect, Thingy
connect("mongodb://localhost/acI")

def load_genomes_to_db():
    im_dict = {
    "IMCC19121" : "IMA7",
    "IMCC26103" : "IMA4",
    "IMCC25003" : "IMA1",
    "IMCC26077" : "IMC1"
    }
    im_ncbi = {
    "IMCC19121" : "NZ_CP015605",
    "IMCC26103" : "NZ_CP015604",
    "IMCC25003" : "NZ_CP015603",
    "IMCC26077" : "NZ_CP015606"
    }
    trans_name = lambda name : "CP" + name.split(".")[0][-2:] if "CP" in name else im_dict[name.split(".")[0]]
    get_ncbi = lambda name : "NZ_" + name.split(".")[0] if "CP" in name else im_ncbi[name.split(".")[0]]
    genome_path = "/home/moritzbuck/people/S001_acI_horizontals/000_data/acI_genomes_annotation/"
    genome_files = [f for f in os.listdir(genome_path) if "gbff" in f]
    for g in genome_files :
        Genome.FromRawFile(file_name = genome_path + g, name=trans_name(g), whatever_id={'ncbi' : get_ncbi(g), 'short' : trans_name(g), 'full' : g}, gtdb_taxon_string=None)
