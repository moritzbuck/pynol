from mongo_thingy import Thingy
from bson.objectid import ObjectId
from pynol.common.COGs.COG import COG
from tqdm import tqdm
from pynol.common.sequence.CDS import CDS
from pynol.common.genome.Genome import Genome
from bson.objectid import ObjectId
import subprocess
import igraph

class COGing( Thingy ):
    """docstring for COG."""

    def __getitem__(self, key):
        if type(key) == int :
            cog_id = self._cogs[key]
            return COG.find_one(ObjectId(cog_id))
        else :
            print("Not implemented yet for non numeric")

    @property
    def cogs(self):
        cogs = list(COG.find({'_id' : {'$in' : self._cogs}}))
        return cogs

    @cogs.setter
    def cogs(self, cogs):
        self._cogs = [f.id for f in cogs]

    @classmethod
    def FromFile(cls, name, file_name, method):
        obj = cls()
        if method == "silix":
            obj.cogs = cls.silix_loader(file_name)
        obj.name = name
        obj.save()
        return obj

    @classmethod
    def FromDict(cls, name, dictionary):
        obj = cls()
        cogs = [COG.FromList(k,v) for k,v in tqdm(dictionary.items())]
        obj.cogs = cogs
        obj.name
        obj.save()
        return obj



    @classmethod
    def silix_loader(cls, file):
        with open(file) as handle:
            pairs = [(l.split()[0], l.strip().split()[1]) for l in handle ]
        cogs = {p[0] : [] for p in pairs}
        for k,v in tqdm(pairs):
            cogs[k] += [ v ]

        cogs_out = []
        for k,v in tqdm(cogs.items()):
            founds = CDS.find({'_id' : {'$in' : [ObjectId(vv) for vv in v]}})
            cc =  COG.FromList(k,list(founds))
            cogs_out += [cc]
        return cogs_out


    @classmethod
    def duck_silix_loader(cls, file):
        with open(file) as handle:
            pairs = [(l.split()[0], l.strip().split()[1]) for l in handle ]
        cogs = {p[0] : [] for p in pairs}
        for k,v in tqdm(pairs):
            cogs[k] += [ v ]

        cogs_clean = {}
        for k,v in tqdm(cogs.items()):
            founds = CDS.find({'_id' : {'$in' : [ObjectId(vv) for vv in v]}})
            cogs_clean[k] = list(founds)
        cogs = [COG.FromList(k,v) for k,v in tqdm(cogs_clean.items())]
        return cogs

    @property
    def cog2genome(self):
        if not self._cog2genome:
            print("pre-computing cog2genome")
            self._cog2genome = {c.id : list(Genome.find({'_id' : {'$in' : list(set(c.genomes))}})) for c in tqdm(self.cogs)}
        return self._cog2genome


    def precluster(self, file, simi = 99):
        if Genome.find_one().species_cluster:
            g2clust = {g.id : g.species_cluster for g in Genome.find()}
            genome_clusters = [[] for x in range(max(g2clust.values())+1)]
            for g in Genome.find() :
                genome_clusters[g.species_cluster] += [g]
        else :
            with open(file) as handle:
                hanis = [h.split() for h in handle.readlines()]
            hanis = { (h[0][:-4], h[1][:-4]) : (float(h[2]), float(h[3]), float(h[4]))  for h in hanis}


            species_graph = igraph.Graph()
            vertexDeict = { v : i for i,v in enumerate(set([x for k in hanis.keys() for x in k]))}
            species_graph.add_vertices(len(vertexDeict))
            species_graph.add_edges([(vertexDeict[k[0]], vertexDeict[k[1]]) for k, v in hanis.items() if float(v[0]) > simi and k[0] != k[1] and float(v[1]) > (float(v[2])*3/4) ])
            vertex2genome = { v : Genome.find_one(ObjectId(k)) for k,v in vertexDeict.items()}
            genome_clusters = [[vertex2genome[cc] for cc in c ] for c in species_graph.components(mode=igraph.WEAK)]

        get_cogs = lambda genome_list : [set([c.id for c in COG.find({'_genomes' : g.id})])  for g in tqdm(genome_list)]
        self.cluster2cog = [list(set.union(*get_cogs(v))) for v in tqdm(genome_clusters)]
        for i, ll in enumerate(genome_clusters):
            for g in ll:
                g.species_cluster = i
                g.save()
        return self.cluster2cog

    def coocurrences(self, cutoff = 0 , genome_set = None, precluster = None):

        dist = lambda v1, v2 : len(v1.intersection(v2))/min(len(v1),len(v2))

        if precluster:
            if not self.cluster2cog:
                tt = self.precluster(precluster)
            all_cogs = set.union(*[set(l) for l in self.cluster2cog])
            c2g_index = { cog : [] for cog in all_cogs}
            for species, cogs in tqdm(enumerate(self.cluster2cog)):
                for cog in cogs :
                    c2g_index[cog] += [species]
            c2g_index = { cog : set(species) for cog, species in c2g_index.items()}
        else :
            c2g_index = {c : set([gg.id for gg in g]) for c,g in self.cog2genome.items()}

        if genome_set :
            genome_set_index = set([g.species_cluster for g in genome_set])
            c2g_index = {k : genome_set_index.intersection(v) for k,v in c2g_index.items()}

        c2g_index = {k : v for k,v in c2g_index.items() if len(v) > cutoff}

        dists = { (k1,k2) : dist(c1,c2) for k1,c1 in tqdm(c2g_index.items()) for k2,c2 in tqdm(c2g_index.items())}
        return dists

    def __iter__(self):
        return COGIterator(self)

    def compute_core(self, cutoff = 0, genome_set = None, name = None, precluster = None):
        rscript = """
            library(pheatmap)
            library(data.table)

            name = "{taxon}"
            matrix = dcast(fread("{taxon}_cooc.csv"), V1~V2)
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
            pheatmap(matrix, cluster_rows=hcl, cluster_cols=hcl, annotation_row=core.dat, main = name, show_colnames=FALSE , show_rownames=FALSE, filename="{taxon}.core.pdf")
            write.table(temp2,"{taxon}.core.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
            """

        dat = self.coocurrences(cutoff, genome_set, precluster = precluster)
        with open("{taxon}_cooc.csv".format(taxon = str(name)), "w") as handle:
            handle.writelines([ str(k[0]) + "," + str(k[1]) + "," +  str(v) + "\n"   for k,v in dat.items()])

        with open("{taxon}.R".format(taxon = str(name)), "w") as handle:
            handle.writelines(rscript.format(taxon = str(name)))
        subprocess.call("Rscript {taxon}.R".format(taxon = str(name)), shell=True)

        with open("{taxon}.core.csv".format(taxon = str(name))) as handle:
            core = set([ObjectId(l[:-1]) for l in handle])

        return core
#        return None

class COGIterator():
    def __init__(self, coging):
        self.i = 0
        self.n = len(coging._cogs)
        self.coging = coging

    def __iter__(self):
        return self

    def __next__(self):
        if self.i < self.n:
            i = self.i
            self.i += 1
            return self.coging[i]
        else:
            raise StopIteration()
