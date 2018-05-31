from pynol.common.genome.sources.Source import Source
from pynol.common.sequence.Genomic import Genomic
from Bio import SeqIO
from Bio.Alphabet import DNAAlphabet
import re
from os.path import splitext

class FromFile(Source):
    """docstring for FromFile."""
    @classmethod
    def FromFile(cls,file_name):
        obj = cls()
        obj.file_name = file_name
        obj.save()
        return obj

    @property
    def file_name(self):
        return self._file_name

    @file_name.setter
    def file_name(self, file_name):
        self._file_name = file_name

    def get_seqio_genome_parser(self) :
        file_name,extension = splitext(self.file_name)
        if extension == ".gbff" :
            typ = "genbank"
        else :
            typ = "fasta"
        return SeqIO.parse(self.file_name, typ, alphabet = DNAAlphabet())
