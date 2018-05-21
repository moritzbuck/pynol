# File: Genomic.py

from pynol.common.sequence.Feature import Feature

class RNA( Feature ):

    @classmethod
    def fromGFFline(cls, line):
        assert "RNA" in line[2] , "This is not an RNA"
        obj = super(RNA,cls).fromGFFline(line)
        obj.save()
        return obj
