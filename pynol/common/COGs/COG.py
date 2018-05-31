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

    @feature_list.setter
    def feature_list(self, features):
        self._feature_list = [f.id for f in features]

    @classmethod
    def FromList(cls, name, feature_list):
        obj = cls()
        assert all([f.id != None for f in feature_list]), "not all Feats are saved in the db"
        obj.feature_list = feature_list
        obj.save()
        obj.name = name
        return obj
