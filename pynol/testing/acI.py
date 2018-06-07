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

connect("mongodb://localhost/acI")

cogs = COGing.find_one()
genomes = list(Genome.find())

def make_line(feature, offset, thresh = 10**-3):
    name = feature.pretty_id.split(":")[0]+"_"+feature.pretty_id.split(":")[-1]
    if offset > 0:
        start = feature.start-offset
        end = feature.end-offset
        strand = feature.strand
    else :
        end = start
    if feature.pfams:
        pfams = ";".join([p.get('query_accession').split(".")[0] for p in feature.pfams if p.get('i-Evalue') < thresh])
    else :
        pfams = ""
    cog = COG.find_one({'_feature_list' : feature.id})
    data = [name, str(start),str(end),str(strand), pfams, cog.name if cog else "", cog.group if cog and cog.group else "", feature.type]
    return ",".join(data)


home = os.environ['HOME']
wei_path = pjoin(home, "people/S001_acI_horizontals/")
data_root = pjoin(wei_path,"000_data/")
temp_root = pjoin(wei_path,"100_temp/")
temp_hmm = pjoin(temp_root,"110_pfams/")
genome_path = pjoin(data_root, "acI_genomes_annotation/")
compliant_path = pjoin(data_root,"compliantFasta/")
group_folder = pjoin(data_root, "cluster_groups")
cluster_file = pjoin(data_root, "groups_with_singletons.txt")
genome_files = [f for f in os.listdir(genome_path) if "gbff" in f]
gi_path = pjoin(wei_path, "200_GI_data")

gi1_island = {g.name : [] for g in genomes}
for c in cogs:
    if c.group:
        if "GI1" in c.group:
            gi1_island[c.group.split(".")[0]] += [c]

for k, v in gi1_island.items():
    feats = sum([vv.feature_list for vv in v ],[])
    feats = [ f for f in feats if Genome.find_one(f.genome).name == k]
    gi1_island[k] = feats

tt = {k : sum([ [vv.start,vv.end] for vv in v],[])  for k,v in gi1_island.items()}
flanks = 40000
genome_ranges = {
                'CP68' : (1138108,1238108),
                'CP79' : (0,100000),
                'IMC1' : (0000,2100000)
                }

gi1_ranges = {k : (min(v)-flanks, max(v)+flanks) if len(v) > 0 else genome_ranges[k]  for k,v in tt.items() }

gi1_features = {}
range_seq_records = []
for k, v in gi1_ranges.items():
    print("Getting all feature in the ranges of GI1 of {genom} within +- {kb}kb".format(genom = k, kb = flanks))
    genome = Genome.find_one({ '_name' : k})
    vv = list(v)
    vv[0] = v[0]#-flanks
    vv[1] = v[1]#+flanks
    gi1_features[k] = Feature.search(genome, vv)


for k, v in gi1_features.items():
    gi1_features[k] = sorted(gi1_features[k],key = lambda bla : bla.start)

maxs ={k : max([ff.start for ff in f]) for k,f in gi1_features.items()}
mins ={k : min([ff.start for ff in f]) for k,f in gi1_features.items()}

for k in gi1_ranges.keys():
    genome = Genome.find_one({ '_name' : k})
    range_seq_records += [SeqRecord(seq = genome.contigs[0].sequence[mins[k]:maxs[k]], id = "GI1_" + genome.name, description= "")]

seq_file = pjoin(temp_root, "120_all_v_all", "GI1_ranges_{kb}kb.fasta".format(kb = flanks/1000))
with open(seq_file, "w") as handle:
    SeqIO.write(range_seq_records,handle,"fasta")

print("blasting")
os.system("makeblastdb -dbtype nucl -out {file} -in {file}".format(file = seq_file))
os.system("tblastx -db {file} -query {file} -outfmt 6 -evalue 0.001 > GI1_v_GI1.tblastx".format(file = seq_file))

for k,v in gi1_features.items():
    file_path = pjoin(gi_path, "210_GI1", k + "_GI1.csv")
    with open(file_path,"w") as handle:
        handle.writelines(["name,start,end,strand,pfams,cog,group,type\n"])
        handle.writelines([make_line(vv, mins[k] ) + "\n" for vv in v])


def parse_pfam_file(file):
    domtblout_head = ["target_name" , "target_accession" , "tlen" , "query_name" , "query_accession" , "qlen" , "E-value","score-sequence" , "bias-sequence" , "nb" , "of" , "c-Evalue" , "i-Evalue" , "score-domain" , "bias-domain" , "hmm_from" , "hmm_to" , "ali_from" , "ali_to" , "env_from" , "env_to" , "acc" , "description_of_target"]
    data = read_csv(file, delim_whitespace=True, comment="#", names = domtblout_head)
    genes = set([dd[1]['target_name'] for dd in data.iterrows()])
    pfam_dict = {ObjectId(g) : [] for g in genes}
    for dd in data.iterrows():
        pfam_dict[ObjectId(dd[1]['target_name'])] += [ dd[1].to_dict() ]
    return pfam_dict

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


## Load genomes
    for g in genome_files :
        print("processing {genome}".format(genome = g))
        Genome.FromRawFile(file_name = genome_path + g, name=trans_name(g), whatever_id={'ncbi' : get_ncbi(g), 'short' : trans_name(g), 'full' : g}, gtdb_taxon_string=None)
        Feature.parse_gff(genome_path + g)

## Load COGs
    compliants = {s.id : s.seq  for g in Genome.find() for s in SeqIO.parse(compliant_path + g.name + ".fasta", "fasta")}
    with open(cluster_file) as handle :
        clusters_lines = { l.split(": ")[0] : l[:-1].split(": ")[1].split() for l in handle}
    clust_sets = {k : [vv for vv in v if vv in compliants] for k,v in clusters_lines.items()}
    clust_sets = { k : v for k,v in clust_sets.items() if len(v) > 0}
    cds_dict = {str(c.protein[:-1]) : c.id  for c in CDS.find()}

    soft_keys = {b : [ k for k in cds_dict.keys() if str(compliants[b])[1:-1] in k ]  for b in  tqdm(sum(list(clust_sets.values()),[]))}
    soft_keys  = { k : v for k,v in soft_keys.items() if len(v) > 0}

    cleaner = {k : [ cds_dict[vv] for vv in v if  Genome.find_one(CDS.find_one(cds_dict[vv]).genome)._name == k.split("|")[0]] for k,v in soft_keys.items()}
    clust_sets = { k : sum([ cleaner.get(vv) for vv in v if cleaner.get(vv) ],[]) for k,v in clust_sets.items() }
    clust_sets = { k : [CDS.find_one(vv) for vv in v] for k,v in clust_sets.items()}
    weis = COGing.FromDict(clust_sets)

    COG_groups = {}
    for p in tqdm(os.listdir(group_folder)):
        with open(pjoin(group_folder, p)) as handle:
            COG_groups[p[:-7]] = [COG.find_one({ 'name' : h.split(": ")[0]}).id for h in handle.readlines()]

    for k, v in COG_groups.items():
        for vv in v:
            cog = COG.find_one(vv)
            cog.group = k
            cog.save()
## run hmms runs

    for g in tqdm(Genome.find()):
        pfam_db = pjoin(home , "data/Pfam-A.hmm")
        hmm_file = pjoin(temp_hmm, g.name + ".hmm")
        prot_file = pjoin(temp_hmm, g.name + ".faa")

        g.write_fasta(prot_file)
        exe_string = "hmmsearch --cpu {cpus} --domtblout {out} {db} {seqs} > /dev/null".format(cpus = 3, out = hmm_file, db = pfam_db, seqs = prot_file)
        os.system(exe_string)

        output = parse_pfam_file(hmm_file)

        for cds, hits in tqdm(output.items()):
            cds_obj = CDS.find_one(cds)
            cds_obj.pfams = hits
            cds_obj.save()

## all_v_all blast

    for g in tqdm(Genome.find()):
        prot_file = pjoin(temp_hmm, g.name + ".fna")
        g.write_fasta(prot_file, pretty = True)
