# File: Source.py
from pynol.common.Record import Record
import re
from mongo_thingy import Thingy
from pynol.common.sequence.Sequence import Sequence
from pynol.common.sequence.CDS import CDS
from bson.objectid import ObjectId
from pynol.common.sequence.Sequence import Sequence

class Source( Thingy ):

# Operations
    def get_data(self, seqset):
        contigs = []
        parsed = [s for s in self.get_seqio_sequence_parser()]
        pads = len(str(len(parsed)))

        for i,s in enumerate(parsed):
            c = Sequence()
            c.sequence = s.seq
            c.seqset = seqset.id
            if not c.__getattribute__('other_ids'):
                c.other_ids = {}
            c.other_ids['original'] = s.id
            c.source = self.id
            c.pretty_id = "temporary_id"
            c.save()
            ori_cds = CDS.find_one(ObjectId(c.other_ids['original']))
            if ori_cds :
                c.genome = ori_cds.genome
                c.ori_cds = ori_cds.id
                self.id_setter = 'lambda c,i, pads : str(c.seqseq) + ";" + str(c.genome) + ";" + str(c.id)'
            c.pretty_id = self.id_setter(c, i, pads)
            c.save()
            contigs += [c]

        return contigs
