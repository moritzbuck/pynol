from pynol.common.sequence.sets.sources.Source import Source
from pynol.common.sequence.Sequence import Sequence
from Bio import SeqIO
import re
from os.path import splitext

class FromFile(Source):
    """docstring for FromFile."""
    @classmethod
    def FromFile(cls,file_name, idsetter = None):
        obj = cls()
        obj.file_name = file_name
        obj.save()
        obj.id_setter = idsetter if idsetter else "lambda seq, i, pads: str(seq._id)"
        return obj

    @property
    def file_name(self):
        return self._file_name

    @file_name.setter
    def file_name(self, file_name):
        self._file_name = file_name

    @property
    def id_setter(self):
        return eval(self._id_setter)

    @id_setter.setter
    def id_setter(self, id_setter):
        self._id_setter = id_setter

    def get_seqio_sequence_parser(self) :
        file_name,extension = splitext(self.file_name)
        if extension == ".gbff" :
            typ = "genbank"
        else :
            typ = "fasta"
        return SeqIO.parse(self.file_name, typ)
