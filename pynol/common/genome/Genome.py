# File: Genome.py
import persistent
from pynol.common.sequence.Genomic import Genomic
from pynol.common.genome.sources.Source import Source
from pynol.common.Record import Record

class Genome( Record ):

# Attributes: Instance

    def __init__(self):
        super(Genome, self).__init__()
        self.source = None  # Source
        self.contigs = None  # Genomic[]
        self.checksum = None

    def populate(self):
        self.contigs = self.source.get_data()
