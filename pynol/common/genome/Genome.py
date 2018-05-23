# File: Genome.py
from pynol.common.sequence.Genomic import Genomic
from pynol.common.genome.sources.Source import Source
from pynol.common.genome.sources.FromRefSeq import FromRefSeq
from pynol.common.genome.sources.FromFile import FromFile

import hashlib
from pynol import PYNOL
from mongo_thingy import Thingy
import re
import sys
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

class Genome( Thingy ):

# Attributes: Instance
    def populate(self):
        self.contigs = self.source.get_data(self)

    @classmethod
    def FromUBAFile(cls,file_name, UBA_id, gtdb_taxon_string):
        genome = cls()
        genome.name = UBA_id
        genome.source = FromFile.FromFile(file_name)
        genome.taxonomy = {}
        genome.taxonomy['gtdb'] = gtdb_taxon_string
        genome.other_ids = {}
        genome.other_ids['uba_id'] = UBA_id
        genome.populate()
        genome.save()
        return genome

    @classmethod
    def FromRefSeq(cls, refseq_id, gtdb_taxon_string, ncbi_taxon_str = None):
        genome = cls()
        genome.name = refseq_id
        genome.source = FromRefSeq.FromRefSeq(s[0])
        genome.taxonomy['ncbi'] = ncbi_taxon_str
        genome.taxonomy['gtdb'] = gtdb_taxon_string
        genome.populate()
        genome.save()
        return genome


    @property
    def source(self):
        source_class = getattr(sys.modules[__name__], self._source_class)
        return source_class.find_one({ '_id' : self._source})

    @source.setter
    def source(self, source):
        source.save()
        self._source = source._id
        self._source_class = source.__class__.__name__


    @property
    def contigs(self):
        return list(Genomic.find({'genome' : self.id}))

    @property
    def checksum(self):
        chks = sorted([f.checksum for f in self.contigs])
        sums = "".join(chks).encode('utf-8')
        return hashlib.md5(sums).hexdigest()

    @property
    def proteins(self):
        if not self.contigs[0]:
            self.predict_proteins()
        return sum([c.proteins for c in self.contigs],[])


    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        assert len(value) < 36, "name is too long, use 'long_name' for this"
        assert bool(re.match('^[A-Za-z0-9_.]+$', value)), "please use only alphanumeric characters of underscore"
        self._name = value

    def save(self) :
        self.source.save()
        if not self._id:
            assert not self.__class__.find_one({'_checksum' : self.checksum}), "This genome is already in the db and has a checksum of {check}".format(check = self.checksum)
            assert not self.__class__.find_one({'_id' : self.pretty_id}), "The genome is already in the db and has an id of {check}".format(check = self.pretty_id)
        Thingy.save(self)

    def predict_proteins(self):
        print("Predicting proteins")

    def to_fasta(self, file):
        to_seq_record = lambda c : SeqRecord(id = c.pretty_id, description = "", seq = c.sequence)
        contigs = [to_seq_record(c) for c in self.contigs]
        with open(file, "w") as handle:
            SeqIO.write(contigs, handle, "fasta")

    def proteom(self, file, pretty = False):
        if pretty:
            to_seq_record = lambda c : SeqRecord(id = c.pretty_id, description = "", seq = c.sequence)
        else :
            to_seq_record = lambda c : SeqRecord(id = c.id, description = "", seq = c.sequence)

        cdss = CDS.find({'genome' : self.id})
        with open(file, "w") as handle:
            SeqIO.write(contigs, handle, "fasta")

    @staticmethod
    def clean(g):
        for gg in Genomic.find({'_pretty_id' : {'$regex' : '{gen}.*'.format(gen = g.name)}}):
            Genomic.delete(gg)
        Genome.delete(g)
