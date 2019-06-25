from mongo_thingy import Thingy
from pynol.common.sequence.Feature import Feature
from bson.objectid import ObjectId

class COG( Thingy ):
    """docstring for COG."""


    @property
    def feature_list(self):
        feats = list(Feature.find({'_id' : {'$in' : self._feature_list}}))
        for cls in Feature.__subclasses__():
            feats += list(cls.find({'_id' : {'$in' : self._feature_list}}))
        return feats

    @property
    def genomes(self):
        if not self._genomes :
            self._genomes = list(set([f.genome for f in self.feature_list]))
            self.save()
        return self._genomes


    @feature_list.setter
    def feature_list(self, features):
        self._feature_list = [f.id for f in features]

    @classmethod
    def FromList(cls, name, feature_list):
        obj = cls()
        assert all([f.id != None for f in feature_list]), "not all Feats are saved in the db"
        obj.feature_list = feature_list
        obj.name = name
        obj.save()
        return obj

    def get_function(self):
        feats = self.feature_list
        data = [f.more for f in feats]
        product = [d.get('product')for d in data]
        hyp_count = product.count('hypothetical protein')
        product = [p for p in product if p != 'hypothetical protein']
        if len(product) > 0:
            best_product = max({s : product.count(s) for s in set(product)}.items(), key = lambda x : x[1])
        else :
            return ('hypothetical protein', 1.0, 1.0)
        product_fraction = best_product[1]/len(feats)
        hyp_fraction = hyp_count/len(feats)
        return (best_product[0] , product_fraction, hyp_fraction)


    def ec_number(self):
        feats = self.feature_list
        data = [f.more for f in feats]
        product = [d.get('eC_number')for d in data]
        hyp_count = product.count(None)
        product = [p for p in product if p != None]
        if len(product) > 0:
            best_product = max({s : product.count(s) for s in set(product)}.items(), key = lambda x : x[1])
        else :
            return ('hypothetical protein', 1.0, 1.0)
        product_fraction = best_product[1]/len(feats)
        hyp_fraction = hyp_count/len(feats)
        return (best_product[0] , product_fraction, hyp_fraction)
