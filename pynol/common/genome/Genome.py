# File: Genome.py
import persistent
from pynol.common.sequence.Genomic import Genomic
from pynol.common.genome.Source import Source

class Genome( persistent.Persistent ):

# Attributes: Instance

    def __init__(self):
        self.source = None  # Source
        self.contigs = None  # Genomic[]
        self.checksum = None 
# Operations
