from mongo_thingy import Thingy
from bson.objectid import ObjectId
from pynol.common.COGs.COG import COG
from tqdm import tqdm
from pynol.common.sequence.CDS import CDS

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
        if method == "duck":
            obj.cogs = cls.duck_silix_loader(file_name)
        obj.name
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
        for k,v in pairs:
            cogs[k] += [v]

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

    def coocurrences(self, cutoff = 0 , genome_set = None):
        dist = lambda v1, v2 : len(v1.intersection(v2))/min(len(v1),len(v2))
        tt = self.cogs
        tt = {c.id : set(c.genomes) for c in tqdm(tt)}
        if genome_set :
            tt = {k : v.intersection(genome_set) for k,v in tt.items()}
        tt = {k : v for k,v in tt.items() if len(v) > cutoff}

        dists = { (k1,k2) : dist(c1,c2) for k1,c1 in tqdm(tt.items()) for k2,c2 in tqdm(tt.items())}
        return dists

    def __iter__(self):
        return COGIterator(self)

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
