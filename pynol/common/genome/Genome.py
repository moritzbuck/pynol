# File: Genome.py
from pynol.common.sequence.Genomic import Genomic
from pynol.common.genome.sources.Source import Source
from pynol.common.genome.sources.FromRefSeq import FromRefSeq
from pynol.common.genome.sources.FromFile import FromFile
from pynol.common.sequence.CDS import CDS
from pynol.common.COGs.COG import COG
from bson.objectid import ObjectId


import hashlib
from pynol import PYNOL
from mongo_thingy import Thingy
import re
import sys
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

class Genome( Thingy ):

# Attributes: Instance

    def __len__(self):
        return self.length

    def populate(self):
        self.source.get_data(self)

    @classmethod
    def FromUBAFile(cls,file_name, UBA_id, gtdb_taxon_string):
        genome = cls()
        genome.name = UBA_id
        genome.source = FromFile.FromFile(file_name)
        genome.taxonomy = {}
        genome.taxonomy['gtdb'] = gtdb_taxon_string
        genome.other_ids = {}
        genome.other_ids['uba_id'] = UBA_id
        genome.save()
        genome.populate()
        return genome

    @classmethod
    def FromRefSeq(cls, refseq_id, gtdb_taxon_string, ncbi_taxon_str = None):
        genome = cls()
        genome.name = refseq_id
        genome.source = FromRefSeq.FromRefSeq(refseq_id )
        genome.taxonomy = {}
        genome.taxonomy['ncbi'] = ncbi_taxon_str
        genome.taxonomy['gtdb'] = gtdb_taxon_string
        genome.save()
        genome.populate()
        return genome

    @classmethod
    def FromRawFile(cls,name, file_name, whatever_id, gtdb_taxon_string, additional_info = {}):
        genome = cls()
        genome.name = name
        genome.source = FromFile.FromFile(file_name)
        genome.taxonomy = {}
        genome.other_ids = {}
        assert type(whatever_id) == dict, "whatever_id needs to be a dictionary with id type as key and ids as value"
        genome.other_ids = whatever_id
        genome.save()
        print("Created Genome with {id}".format(id = genome.id))
        genome.populate()
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
    def length(self):
        return sum([len(l) for l in self.contigs])


    @property
    def checksum(self):
        chks = sorted([f.checksum for f in self.contigs])
        sums = "".join(chks).encode('utf-8')
        return hashlib.md5(sums).hexdigest()

    @property
    def proteins(self):
        cdss = list(CDS.find({ 'genome' : self.id }))
        if len(cdss) == 0:
            print('WARNING : {genome} has no proteins, probably genome prediction has not been run'.format(genome = self.name), file = sys.stderr)
        return cdss


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

    def write_fasta(self, file, pretty = False):
        prots = file.split(".")[-1] == "faa"
        if prots :
            seqs = self.proteins
        else :
            seqs = self.contigs

        if pretty:
            to_seq_record = lambda c : SeqRecord(id = c.pretty_id, description = "", seq = c.protein if prots else c.sequence )
        else :
            to_seq_record = lambda c : SeqRecord(id = str(c.id), description = "", seq = c.protein if prots else c.sequence)

        recs = [to_seq_record(c) for c in seqs]
        with open(file, "w") as handle:
            SeqIO.write(recs, handle, "fasta")

    @staticmethod
    def clean(g):
        for gg in Genomic.find({'_pretty_id' : {'$regex' : '{gen}.*'.format(gen = g.name)}}):
            Genomic.delete(gg)
        Genome.delete(g)

    def completness(self, marker_cogs):
        cogs = [COG.find_one(ObjectId(c)) for c in marker_cogs]
        sc_cogs = [c for c in cogs if len(c.genomes) == len(c._feature_list)]
        hits = [self.id in c.genomes for c in sc_cogs]
        return sum(hits)/len(sc_cogs)
