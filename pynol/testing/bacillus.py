import os
from os.path import join as pjoin
from pynol.common.sequence.Genomic import Genomic
from pynol.common.genome.Genome import Genome
from pynol.common.taxonomy.Taxonomy import Taxonomy
from pynol.common.sequence.Feature import Feature
from tqdm import tqdm
from pynol.common.sequence.CDS import CDS
from pynol.common.sequence.RNA import RNA
from Bio import SeqIO
from pandas import read_csv
from pynol.common.COGs.COG import COG
from pynol.common.COGs.COGing import COGing
from mongo_thingy import connect, Thingy
from bson.objectid import ObjectId
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import re
from pynol.tools.gene_prediction.Prokka import Prokka

def load_data():
    RefSeqs = [k for k in bacili_tax.keys() if not "UBA" in k]
    UBAs = [k for k in bacili_tax.keys() if "UBA" in k]
    for g in tqdm(RefSeqs) :
        Genome.FromRefSeq(g.replace("RS_","").replace("GB_",""), taxonomy.gtdb[g].get_tax_string(full=True) , ncbi_taxon_str = None)
    for g in tqdm(UBAs) :
        uba_folder =  "/home/moritzbuck/uppmax/proj_folder/uppstore2017149/moritz/temp/release80/UBAs/"
        Genome.FromUBAFile(pjoin(uba_folder,g + ".fsa"),g, taxonomy.gtdb[g].get_tax_string(full=True))
    Genome.FromRawFile("Pilis_P2", "/home/moritzbuck/temp/p2.out.contigs.fasta", { 'pilis' : 'P2'}, tt)
    bacies = COGing.FromFile("Bacillus","/home/moritzbuck/people/P01_bacili_pynol/200_clusterings/210_silix/bacilluses.silix","silix")

connect("mongodb://localhost/bacillus")
taxonomy = Taxonomy()
taxonomy.from_gtdb_taxonomy_file("/home/moritzbuck/data/bac_taxonomy_r86.tsv")

bacili_tax = {k : v for k,v in taxonomy.gtdb.items() if v.is_child(taxonomy["g__Bacillus"])}
prok = Prokka.find_one()

same_sub = ["GCA_001938995.1" , "GCF_000715205.1" , "GCF_000714495.2" , "GCF_000691485.2" , "GCF_001704975.1" , "GCF_001468935.1" , "GCF_001675655.1" , "GCF_000828345.1" , "GCF_001183525.1" , "GCF_000172815.1" , "GCF_001444515.1" , "GCF_001578205.1" , "GCF_000017885.4" , "GCF_001307105.1" , "GCF_001038765.1" , "GCF_001038905.1" , "GCF_001043695.1" , "GCF_000935315.1" , "GCF_000604385.1" , "GCF_000444805.1" , "GCF_001429445.1" , "GCF_000715185.1" , "GCF_000444215.1" , "GCF_000828425.1" , "GCF_000691165.1" , "GCF_001677975.1" , "GCF_000972825.1" , "GCF_001766395.1" , "GCA_001566315.1" , "GCF_000828395.1" , "GCF_001938665.1" , "GCA_001938705.1" , "GCF_001938685.1" , "GCF_001286985.1" , "GCF_000828375.1" , "GCF_001653905.1" , "GCF_001518755.1" , "GCF_000800825.1" , "GCF_001482475.1" , "GCA_001038525.1" , "GCF_001017485.1" , "GCF_000300535.1" , "GCF_001857925.1" , "GCF_900119345.1" , "GCF_000949525.1" , "GCF_001420655.1" , "GCA_002421265.1" , "GCF_000508145.1" , "GCF_000789425.2" , "GCF_001029865.1" , "GCF_001578165.1" , "GCF_000831655.1" , "GCF_000299555.1" , "GCF_001687085.1" , "GCF_000353245.1" , "GCF_000590455.1" , "GCF_001038845.1" , "GCF_001039015.1" , "GCF_001426125.1" , "GCF_001191605.1" , "GCF_001548215.1" , "GCF_001908475.1" , "GCF_000972685.1" , "GCF_000701305.1" , "GCF_001265125.1" , "GCA_002382215.1" , "GCF_001038915.1" , "GCF_001038775.1" , "GCF_001043535.1" , "GCF_001043795.1" , "GCF_000986655.1" , "GCF_001700735.1" , "GCF_001431145.1" , "GCF_001431785.1" , "GCF_001543165.1" , "GCF_001896335.1" , "GCF_000828455.1" , "GCF_001444735.1" , "GCF_000691145.1" , "GCA_900094985.1" , "GCF_001457015.1" , "GCF_000264255.1" , "GCF_000225935.1" , "GCF_001278705.1"]

all_bacies_cogs = COGing.find_one()
bacies = list(Genome.find())
all_baci_taxos = {taxonomy[t] : [g for g in bacies if taxonomy[g.taxonomy['gtdb']].is_child(taxonomy[t]) ] for t in set([g.taxonomy['gtdb'] for g in bacies])}


for taxa, genomes  in all_baci_taxos.items():
    if len(genomes) > 50:
        print(taxa)
        matrix = all_bacies_cogs.coocurrences(cutoff=len(genomes) * 0.1, genome_set = set([g.id for g in genomes]))
        with open("test_data/" + taxa.get_tax_string(full =True).split(";")[-1] + ".clstrs", "w") as handle:
            handle.writelines(["COG_1,COG_2,overlap\n"])
            handle.writelines(["{g1},{g2},{dist}\n".format(g1=kk[0], g2=kk[1], dist=vv) for kk,vv in  matrix.items() ])


files = [n for n in  os.listdir("test_data/heatmap_and_cores/") if ".cogs" in n]
