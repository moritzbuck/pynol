# File: Source.py
from pynol.common.Record import Record
import re
from mongo_thingy import Thingy
from pynol.common.sequence.Genomic import Genomic

class Source( Thingy ):

# Operations
    def get_data(self, genome):
        contigs = []
        cdss = []

        for s in self.get_seqio_genome_parser():
            c = Genomic()
            c.sequence = s.seq
            c.type = "contig"
            if not c.__getattribute__('other_ids'):
                c.other_ids = {}
            c.other_ids['original'] = s.id
            contigs += [c]
        # for s in self.seqio_proteome_parser():
        #     c = Genomic()
        #     c.sequence = s.seq
        #     c.type = "contig"
        #     if not c.__getattribute__('other_ids'):
        #         c.other_ids = {}
        #     c.other_ids['original'] = s.id
        #     contigs += [c]
        pads = len(str(len(contigs)))
        for i, s in enumerate(contigs):
            s.pretty_id = "{name}:contig:{i}/{len}".format(name = genome.name, i = str(i+1).zfill(pads), len = len(contigs))
            s.save()
        return contigs
