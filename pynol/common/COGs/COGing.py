from mongo_thingy import Thingy

class COG( Thingy ):
    """docstring for COG."""

    def __init__(self):
        self._cogs = []

    @property
    def feature_list(self):
        feats = list(Feature.find({'_id' : {'$in' : self._feature_list}}))
        for cls in Feature.__subclasses__():
            feats += list(cls.find({'_id' : {'$in' : self._feature_list}}))
        return feats

    @cogs.setter
    def cogs(self, features):
        self._feature_list = [f.id for f in features]

    @classmethod
    def FromFile(cls, name, file_name, method):
        obj = cls()

        obj.feature_list = feature_list
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
            founds = CDS.find({'other_ids.prokka' : {'$in' : v}})
            cogs_clean[k] = list(founds)
        cogs = [COG.FromList(k,v) for k,v in tqdm(cogs_clean.items())]
