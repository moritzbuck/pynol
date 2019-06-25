# File: Genomic.py

from pynol.common.sequence.Feature import Feature
from pynol.common.sequence.Genomic import Genomic
from pynol.common.COGs.COG import COG

class CDS( Feature ):

    @classmethod
    def fromGFFline(cls, line):
        assert line[2] == 'CDS' , "This is not a CDS"
        obj = super(CDS,cls).fromGFFline(line)
        obj.save()
        return obj

    @property
    def protein(self):
        return self.sequence.translate()

    def gff_line(self):
        contig = Genomic.find_one(self.genomic).pretty_id
        source = "prokka"
        typ = "CDS"
        start = str(self.start)
        end = str(self.end)
        score = "."
        strand = "."
        phase = "."
        atts = [ "ID=" + self.pretty_id] + [ k + "=" + v for k,v in self.more.items()]
        if self.cog:
            cc = COG.find_one(self.cog)
            atts += [ "COG_name=" + cc.name  ]#, "COG_fct=" + ":".join([str(v) for v in cc.get_function()])  ]
        atts = ";".join(atts)
        return "\t".join([contig, source, typ, start, end, score, strand, phase, atts]) + "\n"
