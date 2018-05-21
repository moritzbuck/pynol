# File: Genomic.py

from pynol.common.sequence.Feature import Feature

class CDS( Feature ):

    @classmethod
    def fromGFFline(cls, line):
        assert line[2] == 'CDS' , "This is not a CDS"
        obj = super(CDS,cls).fromGFFline(line)
        obj.save()
        return obj

    def get_prot(self):
        return self.sequence.translate()
