from pynol.common.genome.sources.Source import Source
from pynol.common.sequence.Genomic import Genomic
from Bio import SeqIO
from Bio.Alphabet import DNAAlphabet
import re

class FromFile(object):
    """docstring for FromFile."""
    def __init__(self, file, name, long_name = None):
        super(FromFile, self).__init__()
        self.file = file
        self.name = name
        self.long_name = long_name

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        assert len(value) < 36, "name is too long, use 'long_name' for this"
        assert bool(re.match('^[A-Za-z0-9_.]+$', value))
        self._name = value

    def get_data(self):
        contigs = []
        with open(self.file) as handle:
            for s in SeqIO.parse(handle, "fasta", alphabet = DNAAlphabet()):
                c = Genomic()
                c.sequence = s.seq
                c.type = "contig"
                c.other_ids['original'] = s.id
                contigs += [c]

        pads = len(str(len(contigs)))
        for i, s in enumerate(contigs):
            c.pretty_id = "{name}:contig:{i}/{len}".format(name = self.name, i = str(i+1).zfill(pads), len = len(contigs))

        return contigs
