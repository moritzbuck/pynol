# File: Genome.py
import persistent
from pynol.common.sequence.Genomic import Genomic
from pynol.common.genome.sources.Source import Source
from pynol.common.Record import Record
import hashlib


class Genome( Record ):

# Attributes: Instance

    def __init__(self):
        super(Genome, self).__init__()
        self.source = None  # Source
        self.contigs = None  # Genomic[]

    def populate(self):
        self.contigs = self.source.get_data()

    def checksum(self):
        chks = sorted([(f.checksum(), len(f.sequence)) for f in self.contigs])
        sums = "".join().encode('utf-8')
        return hashlib.md5(sums).hexdigest()
