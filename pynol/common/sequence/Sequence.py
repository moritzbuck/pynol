# File: Sequence.py
import hashlib
from Bio.Seq import Seq
from mongo_thingy import Thingy
import re

class Sequence( Thingy ):

    @property
    def sequence(self):
        return Seq(self._sequence)

    @sequence.setter
    def sequence(self, sequence):
        assert sequence.__class__ == Seq
        self._sequence = str(sequence)
        self._checksum = self.checksum

    @property
    def checksum(self):
        return hashlib.md5(self._sequence.encode('utf-8')).hexdigest()

    @property
    def pretty_id(self):
        return self._pretty_id

    @pretty_id.setter
    def pretty_id(self, pretty_id):
        assert bool(re.match('^[A-Za-z0-9:;/_.]+$', pretty_id)), "please use only alphanumeric characters of underscore or dot for names and column or semi-column as separator"
        #assert len(pretty_id) < 37
        self._pretty_id = pretty_id

    def save(self) :
        if not self._id:
#            assert not self.__class__.find_one({'_checksum' : self.checksum}), "This exact sequence is already in the db"
            assert not self.__class__.find_one({'_pretty_id' : self.pretty_id}), "The pretty id {tag} is already in the db".format(tag = self.pretty_id)
        Thingy.save(self)


# Operations
