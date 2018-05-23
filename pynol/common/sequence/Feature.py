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
        obj.type = line[2]
        obj.start = start
        obj.end = end
        obj.pretty_id = genomic.pretty_id+":" + obj.type + ":" + (extras['ID'].split("_")[-1] if extras.get('ID') else line[-1].strip().replace("=","_").replace(" ","_"))
        obj.sequence = seq.reverse_complement() if rev_comp else seq
        obj.genomic = genomic.id
        obj.genome = genomic.genome
        obj.other_ids={}
        if extras.get('locus_tag'):
            obj.other_ids['prokka'] = extras['locus_tag']
            del extras['locus_tag']
        if extras.get('ID'):
            del extras['ID']

        obj.more = extras

        obj.save()
        return obj
