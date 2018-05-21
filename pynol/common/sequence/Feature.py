from pynol.common.sequence.Nucleotide import Nucleotide
from pynol.common.sequence.Genomic import Genomic
from Bio.Alphabet import DNAAlphabet

class Feature( Nucleotide ):

    @classmethod
    def fromGFFline(cls, line):
        obj = cls()
        genomic = Genomic.find_one({'_pretty_id' : line[0]})
        start = int(line[3])
        end = int(line[4])
        rev_comp = line[6] == '-'
        seq = genomic.sequence[(start-1):end]
        extras = {f.split("=")[0] : f.split("=")[1]  for f in line[-1].split(";")}
        seq.alphabet = DNAAlphabet()
        obj.pretty_id = genomic.pretty_id+":CDS:" + extras['ID'].split("_")[1]
        obj.sequence = seq.reverse_complement() if rev_comp else seq
        obj.genomic = genomic.id
        obj.other_ids={}
        obj.other_ids['prokka'] = extras['locus_tag']
        obj.type = line[2]

        del extras['ID']
        del extras['locus_tag']

        obj.more = extras

        obj.save()
        return obj
