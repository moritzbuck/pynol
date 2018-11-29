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
from numpy import ceil
import igraph
import subprocess
from mongo_thingy import connect, Thingy
import pandas
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

connect("mongodb://localhost/pateci")

all_pates = list(Genome.find({'analysis' : 'pates'}))
coging = COGing.find_one({'name' : 'patescis'})

taxonomy = Taxonomy()
taxonomy.from_gtdb_taxonomy_file("/home/moritzbuck/data/gtdbtk_taxonomy.tsv")

with open("patescis.gtdbtk") as handle:
    temp = handle.readlines()
temp = [t[:-1].split("\t")[0:2] for t in temp]
genome2tax =  {t[0] : taxonomy.taxa_dict[[tt for tt in t[1].split(";") if not tt.endswith("__") and not tt.startswith("s__") ][-1]]  for t in temp[1:]}

with open("all_v_all.fastani") as handle:
    hanis = [h.split() for h in handle.readlines()]
hanis = { (h[0][:-4], h[1][:-4]) : (float(h[2]), float(h[3]), float(h[4]))  for h in hanis}

extras = set([str(g.id) for g in all_pates]).difference(set(sum([[a[0],a[1]] for a in hanis], [])))

species_graph = igraph.Graph()
vertexDeict = { v : i for i,v in enumerate(set([x for k in hanis.keys() for x in k]))}
species_graph.add_vertices(len(vertexDeict))
species_graph.add_edges([(vertexDeict[k[0]], vertexDeict[k[1]]) for k, v in hanis.items() if float(v[0]) > 95 and k[0] != k[1] and float(v[1]) > (float(v[2])*3/4) ])
vertex2genome = { v : Genome.find_one(ObjectId(k)) for k,v in vertexDeict.items()}
genome_clusters = [[vertex2genome[cc] for cc in c ] for c in species_graph.components(mode=igraph.WEAK)]+ [ [Genome.find_one(ObjectId(g))] for g in extras]


temp = [list(set([genome2tax[str(vv.id)] for vv in v] )) for v in genome_clusters]

specials = [i for i,f in enumerate(genome_clusters) if len(f) > 1 and not taxonomy.check_monophyly([genome2tax[str(ff.id)] for ff in f])]

temp = [ max(t,key = lambda t : t.level) if len(t) > 1 else t[0] for t in temp ]
assert len(specials) == 1
temp[specials[0]] = temp[specials[0]].parent

genome_clust_taxa = temp
uniq_taxa = set([t for t in genome_clust_taxa])
uniq_count = {t : 1 for t in uniq_taxa}
genome_clust_ids = [t for t in genome_clust_taxa]

for i,taxon in enumerate(genome_clust_ids):
    itera = uniq_count[taxon]
    uniq_count[taxon] += 1
    genus = taxon.tail().startswith("g__")
    genome_clust_ids[i] = taxon.tail() + ("_sp_cluster_" if genus else "_new_" + Taxon.rev_leves[taxon.level+1]) + "_" + str(itera)

genome_clusters = { v :k  for k,v in zip(genome_clusters, genome_clust_ids)}
cluster2taxon = { v :k  for k,v in zip(genome_clust_taxa, genome_clust_ids)}
all_taxa = set([taxonomy[t.get_parent(i,full=True)] for t in cluster2taxon.values() for i in range(2,6) if t.get_parent(i,full=True)])

tax2cluster = { tax : [g for g, tt in cluster2taxon.items() if tt.is_child(tax) ] for tax in all_taxa}

big_clades = {k : v for k,v in tax2cluster.items() if len(v) > 20 }

cog2genome = coging.cogs
cog2genome = {c.id : c.genomes for c in tqdm(cog2genome)}
cog2genomes_set = { c    : set(g) for c, g in cog2genome.items()}

genome2cluster = {str(vv.id) : k for k, v in genome_clusters.items() for vv in v}

cog2cluster = {c :[] for c in cog2genome}
genome2cog = {g.id : [] for g in all_pates}
for c, g in tqdm(cog2genomes_set.items()):
    for gg in g:
        genome2cog[gg] +=[c]
        cog2cluster[c] += [genome2cluster[str(gg)]]


genome2cog = {str(c) : set(g) for c, g in genome2cog.items()}
cog2cluster = {str(c) : set(g) for c, g in cog2cluster.items()}
cluster2cog = {c : [ cog for cog, species in cog2cluster.items() if c in species] for c in tqdm(set(cluster2taxon.keys()))}



dist = lambda v1, v2 : len(v1.intersection(v2))/min(len(v1),len(v2))

cluster_v_cluster = {(cls1,cls2) : dist(set(cogs1), set(cogs2)) for cls1,cogs1 in tqdm(cluster2cog.items()) for cls2,cogs2 in cluster2cog.items()}
with open("test_data/cluster_v_cluster.csv", "w") as handle:
        handle.writelines([ str(k[0]) + "," + str(k[1]) + "," + str(v) + "\n" for  k,v in cluster_v_cluster.items()])


for t, ll in tqdm(big_clades.items()):
    cogs = set.union(*[set(cluster2cog[g]) for g in ll])
    cutoff=ceil(len(ll)*0.1)
    cc2cc = { c : set(cog2cluster[c]).intersection(ll) for c in cogs}
    cc2cc = {c : v for c, v in cc2cc.items() if len(v) > cutoff}
    coocs = {(k,l) : dist(v,w) for k, v in tqdm(cc2cc.items()) for l, w in cc2cc.items()}
    with open("test_data/coocs/" + t.tail() + ".csv", "w") as handle:
        handle.writelines([ str(k[0]) + "," + str(k[1]) + "," + str(v) + "\n" for  k,v in coocs.items()])


rscript = """
library(pheatmap)
library(data.table)

name = "{taxon}"
matrix = dcast(fread("test_data/coocs/{taxon}.csv"), V1~V2)
rnames = matrix$V1
matrix = as.data.frame(matrix[, V1 := NULL])
row.names(matrix) = rnames
hcl = hclust(dist(matrix))
clst = data.frame(core = as.factor(cutree(hcl,15)))
mean.overlap = sapply(1:15, function(x) mean(rowMeans(matrix[clst == x,])))
core = which.max(mean.overlap)
core.dat = data.frame(core = (clst == core )+1-1)
temp2 = row.names(core.dat)[core.dat$core == 1]
core.dat$core = as.factor(core.dat$core)
pheatmap(matrix, cluster_rows=hcl, cluster_cols=hcl, annotation_row=core.dat, main = name, show_colnames=FALSE , show_rownames=FALSE, filename="test_data/cores/{taxon}.core.pdf")
write.table(temp2,"test_data/cores/{taxon}.core.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
"""

for t in tqdm(big_clades.keys()):
    with open("temp.R" , "w") as handle:
        handle.writelines(rscript.format(taxon = t.tail()))
    subprocess.call("Rscript temp.R", shell=True)

tax2cores = {}
for t in tqdm(big_clades.keys()):
    with open("test_data/cores/{taxon}.core.csv".format(taxon = t.tail())) as handle:
        tax2cores[t] = set([l[:-1] for l in handle])

all_cores = set.union(*tax2cores.values())

root = min(tax2cores.keys(), key = lambda t : t.level)
hierarch_core = {}
cumul_core = {}
hierarch_core[root] = tax2cores[root]
cumul_core[root] = tax2cores[root]

levels = [i+1 for i in range(root.level, max(tax2cores.keys(), key = lambda t : t.level).level)]
for l in levels:
    for t, core in tax2cores.items():
        if t.level == l:
            hierarch_core[t] = core.difference(cumul_core[t.parent])
            cumul_core[t] = core.union(cumul_core[t.parent])

pres_abs_mat = {cog : {cluster2taxon[clust].get_tax_string(full = True) + ";" + clust : int(cog in genes) for clust, genes in cluster2cog.items()} for cog in tqdm(all_cores)}

pandas.DataFrame.from_dict(pres_abs_mat, orient="index").to_csv("test_data/core_matrix.csv")
pandas.DataFrame({ Taxon.rev_leves[i] : {  cluster2taxon[k].get_tax_string(full = True) + ";" + k: t.get_parent(i) for k, t in cluster2taxon.items()} for i
in range(2,7)}).to_csv("test_data/cluster_metadata.csv")

core_mat = {cog : {clust.get_tax_string(full=True) : int(cog in genes) for clust, genes in cumul_core.items()} for cog in tqdm(all_cores)}
pandas.DataFrame.from_dict(core_mat,orient="index").to_csv("test_data/cog_md.csv")


all_core_cogs = [COG.find_one(ObjectId(c)) for c in set.union(*tax2cores.values())]

cog2function = { cog.id : cog.get_function()[0] for cog in tqdm(all_core_cogs)}
function2cog = {f : [c for c, ff in cog2function.items() if ff == f] for f in set(cog2function.values())}
function2cog['hypothetical protein'] = []


for taxon, core in tqdm(hierarch_core.items()):
    data = {}
    cluster_set = tax2cluster[taxon]
    genome_set=[g for g in sum([genome_clusters[c]  for c in tax2cluster[taxon]],[])]
    cogs = [COG.find_one(ObjectId(c)) for c in core]
    temp_dict = {g.id : [] for g in genome_set}
    for cog in tqdm(cogs):
        funct = cog.get_function()
        feats = cog.feature_list
        count = 0
        for f in feats:
            if temp_dict.get(f.genome) != None:
                count += 1
                temp_dict[f.genome] += [f]
        found_genome = [Genome.find_one(g).name for g in cog.genomes if genome2tax[str(g)].is_child(taxon)]
        within = cog2cluster[str(cog.id)].intersection(cluster_set)
        without = cog2cluster[str(cog.id)].difference(cluster_set)
        found_fraction = len(within)/len(cluster_set)
        found_outside = len(without)/len(genome_clusters)
        enrichment = found_fraction/found_outside if found_outside != 0 else float("inf")
        data[cog.id] = {
                        'name' : cog.name,
                        'function' : funct[0],
                        'fct_fract' : funct[1],
                        'hypo_fract' : funct[2],
                        'found_fraction' : found_fraction,
                        'found_outside' : found_outside,
                        'enrichment' :  enrichment,
                        'genome_count' : len(found_genome),
                        'copy_per_genome' : count/len(genome_set),
                        'functional_similar_cogs': ";".join([COG.find_one(c).name for c in function2cog[funct[0]] if c != cog.id])
                        }


    cog_clustering = make_proximity_graph(temp_dict)
    for cog, group in cog_clustering.items():
        data[cog]['operon'] = group
    pandas.DataFrame.from_dict(data, orient='index').sort_values(by="operon").to_csv(pjoin('test_data/cores/', taxon.tail() + "_core_data.csv"))

for cog in tqdm(all_core_cogs):
    SeqIO.write([SeqRecord(seq = c.protein, name = "", id =
"{tax}:{pid}:{fct}".format(tax = genome2tax[str(c.genome)].get_tax_string().replace(';','.'), pid = c.pretty_id, fct = c.more['product']), description = "") for c in cog.feature_list], pjoin('test_data/cogs/', "COG_" + cog.name + ".fasta") , "fasta")




all_cooc_cogs = set.union(*[set.union(*s) for s in grafs.values()])

for g in tqdm(all_pates):
    g.write_fasta("/home/moritzbuck/temp/all_patesc/{id}.faa".format(id=g.id))

for g in tqdm(all_pates):
    g.write_fasta("/home/moritzbuck/temp/all_patesc/{id}.fna".format(id=g.id))



def make_proximity_graph(temp_dict, dist_cutoff = 10000, count_cutoff = 0.3):
    all_cogs = set([f.cog for f in sum(temp_dict.values(),[])])

    all_edges = {(f,g) : 0 for f in all_cogs for g in all_cogs}
    for k, v in temp_dict.items():
        for t in v:
            for z in v:
                if t.genomic == z.genomic and abs(t.start-z.start) < dist_cutoff and t != z:
                    all_edges[(t.cog, z.cog)] += 1

    abundant_edges = {k : v for k, v in all_edges.items() if v > len(temp_dict)*count_cutoff}

    prox_graph = igraph.Graph()
    vertexDeict = { v : i for i,v in enumerate(set(sum([list(x) for x in abundant_edges.keys()], [])))}
    prox_graph.add_vertices(len(vertexDeict))
    prox_graph.add_edges([(vertexDeict[x[0]], vertexDeict[x[1]]) for x in abundant_edges.keys()])
    vertex2cog = { v : k for k,v in vertexDeict.items()}
    cog_grafs = [[vertex2cog[cc] for cc in c ] for c in prox_graph.components(mode=igraph.WEAK)]

    return {cc : str(i) if len(c) > 1 else 'None' for i,c in enumerate(cog_grafs) for cc in c }
