from ete3 import Tree, NodeStyle, TreeStyle
from ete3 import COLOR_SCHEMES
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
from scipy.stats import hypergeom
from pandas import DataFrame
import statsmodels.stats.multitest as multi
from Bio.SeqRecord import SeqRecord
from pynol.tools.similarity_search.Blast import Blast
from Bio import Entrez
import re

connect("mongodb://localhost/acI")

cogs = COGing.find_one()
genomes = list(Genome.find())

with open("/home/moritzbuck/data/pfamA.txt") as handle:
    pfam2name={re.split(r'\t',l)[0] : re.split(r'\t',l)[3] for l in handle.readlines()}


with open("/home/moritzbuck/data/pfamA.txt") as handle:
    pfam2name={re.split(r'\t',l)[0] : re.split(r'\t',l)[3] for l in handle.readlines()}


def make_line(feature, offset, thresh = 10**-3):
    name = feature.pretty_id.split(":")[0]+"_"+feature.pretty_id.split(":")[-1]
    product = feature.more.get('product')
    if offset > 0:
        start = feature.start-offset
        end = feature.end-offset
        strand = feature.strand
    else :
        end = start
    if feature.pfams:
        pfams = ";".join([p.get('query_accession').split(".")[0] for p in feature.pfams if p.get('i-Evalue') < thresh])
        pfam_names = ";".join([pfam2name[p.get('query_accession').split(".")[0]] for p in feature.pfams if p.get('i-Evalue') < thresh])
    else :
        pfams = ""
        pfam_names = ""
    cog = COG.find_one({'_feature_list' : feature.id})
    data = [name, product, str(start),str(end),str(strand), pfams, pfam_names , cog.name if cog else "", ";".join(cog.group) if cog and cog.group else "", feature.type]
    return "\t".join(data)


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
ancestor_path = pjoin(gi_path, "220_ancestor_gains")

[os.makedirs( pjoin(ancestor_path , c.name )) for c in cogs if c.group and "ancestor_gains" in c.group if not os.path.exists( pjoin(ancestor_path , c.name ))]
[ SeqIO.write([ SeqRecord(id = c.name + ";" + Genome.find_one(cc.genome).name , seq =cc.protein, description = "") for cc in c.feature_list ], pjoin(ancestor_path , c.name , c.name + ".faa"), "fasta") for c in cogs if c.group and "ancestor_gains" in c.group]

os.system("cat {path}/*/*.faa > {path}/all_gained.faa".format(path=ancestor_path))
#os.system("blastp  -num_threads 20  -query all_gained.faa  -db /sw/data/uppnex/blast_databases/nr -outfmt 10  >  {path}/all_gained_vs_nr.blastp".format(path = ancestor_path))

raw_all_blast = Blast.parse_blast_file("/home/moritzbuck/people/S001_acI_horizontals/200_GI_data/220_ancestor_gains/all_gained_vs_nr.blastp")
raw_all_blast['cog'] = [c.split(";")[0] for c in raw_all_blast['query']]

cutoffs = dict(raw_all_blast.groupby('cog').apply(lambda df : {ii: len(set(df.loc[df['evalue'] < 10**-ii].subject)) for ii in list(range(3,10))+[20,30,40,50] } ))
min_seq = 250
min_cutoff = lambda dd, v : min({ k : vv for k, vv in dd.items() if vv > min_seq }.items(), key=lambda x : x[1])
cutoffs = {cog : (min_cutoff(v, min_seq) if max(v.values()) > min_seq else max(v.items(), key=lambda t : t[1]))[0] for cog, v in cutoffs.items()}
hits = dict(raw_all_blast.groupby('cog').apply(lambda df : set(df.loc[df['evalue'] < 10**-cutoffs[list(df['cog'])[0]]].subject) ))

from Bio import Entrez
Entrez.email = "murumbii@gmail.com"

homos = dict()
for k, ids in tqdm(hits.items()):
    if not homos.get(k):
        with Entrez.efetch(db="protein", rettype="gp", retmode="text", id=list(ids)) as handle:
            homos[k] = [s for s in SeqIO.parse(handle, "genbank")]

with Entrez.efetch(db="protein", rettype="gp", retmode="text", id=bla) as handle:
    hh = [s for s in tqdm(SeqIO.parse(handle, "genbank"))]



[ SeqIO.write([ SeqRecord(id = c.name + ";" + Genome.find_one(cc.genome).name , seq =cc.protein, description = "") for cc in c.feature_list ] + hits[c], pjoin(ancestor_path , c.name , c.name + "_with_hits.faa"), "fasta") for c in cogs if c.group and "ancestor_gains" in c.group]

slurm_script = """#!/bin/bash
#SBATCH -D /home/moritz/temp/trees/
#SBATCH -J trees_{name}
#SBATCH -o /home/moritz/temp/trees/{name}.out
#SBATCH -e /home/moritz/temp/trees/{name}.err
#SBATCH -A snic2017-1-616
#SBATCH -t 10-00:00:00
#SBATCH -n 1
#SBATCH -p core

mkdir {name}
cd {name}
mv ../../{name}_with_hits.faa .
muscle -in {name}_with_hits.faa -out {name}_aligned.faa
fasttree < {name}_aligned.faa > {name}.tree
"""

raxml_slurm_script = """#!/bin/bash
#SBATCH -D /home/moritz/temp/more_trees/
#SBATCH -J trees_{name}
#SBATCH -o /home/moritz/temp/more_trees/rax_{name}.out
#SBATCH -e /home/moritz/temp/more_trees/rax_{name}.err
#SBATCH -A snic2017-1-616
#SBATCH -t 2-00:00:00
#SBATCH -n 4
#SBATCH -p core

module load bioinfo-tools
module load raxml

mkdir -p {name}
cd {name}
mv ../fastas/rax_{name}.faa .
muscle -in rax_{name}.faa -out {name}_aligned.faa
raxmlHPC-PTHREADS-AVX -m PROTGAMMALGF -T 4 -p 42 -s {name}_aligned.faa -n tree -n {name} -f a -x 1 -N autoMR
"""


for c in cogs:
    if c.group and "ancestor_gains" in c.group:
        scr = slurm_script.format(name = c.name)
        with open(pjoin(temp_root, c.name + ".sh"), "w") as handle:
            handle.writelines(scr)

for c in cogs:
    if c.group and "ancestor_gains" in c.group:
        prin("making tree for ", c.name)
        make_tree_fig(pjoin(ancestor_path, c.name, c.name + ".tree"), pjoin(temp_root, c.name, c.name + ".tree"))



gi1_island = {g.name : [] for g in genomes}
for c in cogs:
    if c.group:
        grps = [g for g in c.group if "GI1" in g]
        for g in grps:
            gi1_island[g.split(".")[0]] += [c]

for k, v in gi1_island.items():
    feats = sum([vv.feature_list for vv in v ],[])
    feats = [ f for f in feats if Genome.find_one(f.genome).name == k]
    gi1_island[k] = feats

tt = {k : sum([ [vv.start,vv.end] for vv in v],[])  for k,v in gi1_island.items()}
flanks = 20000
genome_ranges = {
                'CP68' : (1160000,1240000),
                'CP82': (1230030, 1302615),
                'CP79' : (0,50000),
                'IMC1' : (100000,140000)
                }

gi1_ranges = {k : genome_ranges[k] if genome_ranges.get(k) else (min(v)-flanks, max(v)+flanks)  for k,v in tt.items() if len(v) > 0 or genome_ranges.get(k)}

gi1_features = {}
range_seq_records = []
for k, v in gi1_ranges.items():
    print("Getting all feature in the ranges of GI1 of {genom} within +- {kb}kb".format(genom = k, kb = flanks))
    genome = Genome.find_one({ '_name' : k})
    vv = list(v)
    vv[0] = v[0]#-flanks
    vv[1] = v[1]#+flanks
    feats =  Feature.search(genome, vv)
    assert len(feats.values()) == 1
    gi1_features[k] = list(feats.values())[0]

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


full_pfams = {}
for g in gi1_ranges:
    genome = Genome.find_one({ '_name' : k})
    full_pfams[g] = sum([[f['query_accession'] for f in f.pfams if f['E-value'] < 0.001 ]  for f in list(Feature.search(genome).values())[0] if f.pfams],[])

genome2pfams = full_pfams
pfam2genome = {f : [k for k,v in genome2pfams.items() if f in v] for f in full_pfams }

full_pfams = sum(full_pfams.values(),[])
all_pfams = set(full_pfams)
gi1_pfams = sum([[f['query_accession'] for f in f.pfams if f['E-value'] < 0.001 ] for f in sum(gi1_features.values(),[]) if f.pfams],[])

def pfam_hyg(pfam):
        k = gi1_pfams.count(pfam)
        M = len(full_pfams)
        N = len(gi1_pfams)
        n = full_pfams.count(pfam)
        p = hypergeom.sf(k = k, M = M, n = n, N = N)
        ratio = (float(k)/N)/(float(n)/M)
        return p,ratio

pfams_pvals = {p : pfam_hyg(p)[0] for p in all_pfams}
pfams_effect = {p : pfam_hyg(p)[1] for p in all_pfams}

adj_pval = dict(zip(all_pfams,multi.multipletests(list(pfams_pvals.values()), method="fdr_bh")[1]))
sigs = { k.split(".")[0] : { 'name' : pfam2name[k.split(".")[0]], 'pval' : pfams_pvals[k] , 'adj.pval' : v, 'ratio' :pfams_effect[k]} for k,v in sorted(adj_pval.items(), key = lambda x : x[1]) if v < 0.05}
pfam_table = DataFrame.from_dict(sigs, orient = 'index')
pfam_table = pfam_table.sort_values(by = 'adj.pval')
pfam_table['group'] = ["other transferase" if "ransferase" in p else "" for p in pfam_table.name]
pfam_table['group'] = ["methyltransferase" if "Methyl" in p and "transferase" in p else g for p,g in zip(pfam_table.name, pfam_table.group)]
pfam_table['group'] = ["glycosyl transferase" if "Glycosyl" in p and "transferase" in p else g for p,g in zip(pfam_table.name, pfam_table.group)]
pfam_table['group'] = ["replication" if  "DNA" in p or  "elicase"  in p or "opoisomerase" in p else g for p,g in zip(pfam_table.name, pfam_table.group)]
pfam_table.to_csv("/home/moritzbuck/people/S001_acI_horizontals/pfam_table.csv")



print("blasting")
os.system("makeblastdb -dbtype nucl -out {file} -in {file}".format(file = seq_file))
os.system("tblastx -db {file} -query {file} -outfmt 6 -evalue 0.001 > GI1_v_GI1.tblastx".format(file = seq_file))

for k,v in gi1_features.items():
    file_path = pjoin(gi_path, "210_GI1", k + "_GI1.csv")
    with open(file_path,"w") as handle:
        handle.writelines(["name\tproduct\tstart\tend\tstrand\tpfams\tpfamNames\tcog\tgroup\ttype\n"])
        handle.writelines([make_line(vv, mins[k] ) + "\n" for vv in v])

all_gi_cdss = [gg.id for g in gi1_features.values() for gg in g]

links = {cds : sum([ [CDS.find_one(f) for f in COG.find_one(c)._feature_list  if f in all_gi_cdss ]  for c in cds.cogs ],[]) for island in gi1_features.values() for cds in island if cds.cogs}

gi1_mins = { k : min([f.start for f in fs])for k, fs in  gi1_features.items()}
def line_maker(cds):
    nome = Genome.find_one(cds.genome).name
    offset = gi1_mins[nome]
    group = ":".join(sum([COG.find_one(cog).group for cog in cds.cogs if COG.find_one(cog).group],[]))
    return "\t".join([Genome.find_one(cds.genome).name, str(cds.start-offset), str(cds.end-offset), group])

with open("/home/moritzbuck/people/S001_acI_horizontals/gi1_links.txt", "w") as handle:
     handle.writelines([line_maker(k) + "\t" + line_maker(ll) + "\n" for k, l in links.items() for ll in l])



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
            if type(cog.group) == str:
                del cog.group
            if not cog.group:
                cog.group = []
            cog.group += [k]
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


    def make_tree_fig(tree_file, out_name, tax_level = None ):
        with open(tree_file) as handle:
            lines = handle.readlines()

            if len(lines) > 0 :
                ete_tree = Tree(lines[0][:-1].replace(";IM","-IM").replace(";CP","-"))
            else :
                return None
        if tax_level:
            taxa = {}
#            taxa = {xx : [x for x in xx.name.replace(" ","-").split("_")[1:] if len(x) > 0 and not x[0].isdigit() ] for xx in ete_tree.get_leaves()}
            for xx in ete_tree.get_leaves():
                id = xx.name
                taxon = taxas.get(id)
                if taxon :
                    xx.name = ";".join([ id ] + taxon if taxon else [])
                taxa[xx] = taxon[tax_level] if taxon and len(taxon) > tax_level else None


            for leaf in taxa:
                leaf.set_style(NodeStyle())
                if taxa.get(leaf) and cols.get(taxa[leaf]):
                    leaf.img_style["bgcolor"] = cols[taxa[leaf]]
                elif "acI" in leaf.name:
                        leaf.img_style["bgcolor"] = cols['acI']
        else :
            taxa = None
        styl = TreeStyle()
        styl.mode = 'c'
#        styl.arc_start = -180
#        styl.arc_span = 180 #
        print(out_name)
        ete_tree.render(out_name,w=len(ete_tree.get_leaves())*5, tree_style= styl)




cp69tocp77 = {k : [vv for vv in v if "CP77" in vv.pretty_id] for k, v in links.items() if "CP69" in k.pretty_id}
ll = list({ COG.find_one(k.cogs[0]).name : COG.find_one(k.cogs[0]).name for k,v in cp69tocp77.items() if len(v) > 0 }.keys())


[os.makedirs( pjoin(temp_root, "130_other" , c.name )) for c in cogs if c.name in ll and not os.path.exists( pjoin(temp_root, "130_other" , c.name ))]
[ SeqIO.write([ SeqRecord(id = str(i).zfill(2) + "_" + c.name + "_" + Genome.find_one(cc.genome).name , seq =cc.protein, description = "") for i, cc in enumerate(c.feature_list) ], pjoin(temp_root, "130_other" , c.name , "rax_" +  c.name + ".faa"), "fasta") for c in cogs if c.name in ll]

os.system("cat {path}/*/*.faa > {path}/all_gained.faa".format(path=pjoin(temp_root, "130_other")))
#os.system("blastp  -num_threads 20  -query all_gained.faa  -db /sw/data/uppnex/blast_databases/nr -outfmt 10  >  {path}/all_gained_vs_nr.blastp".format(path = ancestor_path))

raw_all_blast = Blast.parse_blast_file("/home/moritzbuck/uppmax/temp/all_gained_vs_nr.blastp")
raw_all_blast['cog'] = [c.split(";")[0] for c in raw_all_blast['query']]

cutoffs = dict(raw_all_blast.groupby('cog').apply(lambda df : {ii: len(set(df.loc[df['evalue'] < 10**-ii].subject)) for ii in list(range(3,10))+[20,30,40,50] } ))
min_seq = 250
min_cutoff = lambda dd, v : min({ k : vv for k, vv in dd.items() if vv > min_seq }.items(), key=lambda x : x[1])
cutoffs = {cog : (min_cutoff(v, min_seq) if max(v.values()) > min_seq else max(v.items(), key=lambda t : t[1]))[0] for cog, v in cutoffs.items()}
hits = dict(raw_all_blast.groupby('cog').apply(lambda df : set(df.loc[df['evalue'] < 10**-cutoffs[list(df['cog'])[0]]].subject) ))

from Bio import Entrez
Entrez.email = "murumbii@gmail.com"

homos = dict()
for k, ids in tqdm(hits.items()):
    if not homos.get(k):
        with Entrez.efetch(db="protein", rettype="gp", retmode="text", id=list(ids)) as handle:
            homos[k] = [s for s in SeqIO.parse(handle, "genbank")]


[ SeqIO.write([ SeqRecord(id = c.name + ";" + Genome.find_one(cc.genome).name , seq =cc.protein, description = "") for cc in c.feature_list ] + homos[c.name], pjoin(temp_root, "130_other" , c.name , c.name + "_with_hits.faa"), "fasta") for c in cogs if c.name in ll]
taxas = {s.id : s.annotations['taxonomy'] for s in sum(homos.values(), [])}

all_phyla = [t[1] for t in taxas.values() if len(t) > 2]
all_phyla_sub = [s for s in set(all_phyla) if all_phyla.count(s)  > 300] + ['acI']



[ SeqIO.write([ SeqRecord(id = str(i).zfill(2) + "_" + c.name + "_" + Genome.find_one(cc.genome).name , seq =cc.protein, description = "") for i, cc in enumerate(c.feature_list) ], pjoin(temp_root, "130_other" , "more_trees", "fastas" , "rax_" +  c.name + ".faa"), "fasta") for c in cogs]

for c in cogs:
#    if c.name in ll:
        scr = raxml_slurm_script.format(name = c.name)
        with open(pjoin(temp_root, "130_other", "more_trees", "scripts",  "rax_" + c.name + ".sh"), "w") as handle:
            handle.writelines(scr)

for c in cogs:
    if c.name in ll and not c.name == 'acIfull00695' :
        print("making tree for ", c.name)
        make_tree_fig(pjoin(temp_root, "130_other", c.name, c.name + ".tree"), pjoin(temp_root, "130_other", c.name, c.name + ".pdf"))
