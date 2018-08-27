from pynol.common.sequence.Nucleotide import Nucleotide
from pynol.common.sequence.Genomic import Genomic
from Bio.Alphabet import DNAAlphabet
from Bio.GenBank import parse
from bson.objectid import ObjectId
import pynol

class Feature( Nucleotide ):
    @classmethod
    def search(cls, genomic, interval = None):
        if type(genomic) == pynol.common.genome.Genome.Genome :
            feats = {}
            for contig in Genomic.find({'genome' : genomic._id}):
                feats[contig] = Feature.search(contig, interval = interval)
            return feats
        else :
            if type(genomic) == ObjectId :
                genomic_query = { 'genomic' : genomic }
            else :
                genomic_query = {'genomic' : genomic.id}
            if interval :
                start_query = { 'start' : {'$gte' : interval[0] } }
                end_query = { 'end' : {'$lte' : interval[1] } }
                query = { '$and' : [ genomic_query , start_query, end_query]}
            else :
                query = genomic_query

            feats = list(cls.find(query))
            for sub_cls in cls.__subclasses__():
                feats += list(sub_cls.find(query))
            return feats

    @classmethod
    def parse_gff(cls, file):
        from pynol.common.sequence.CDS import CDS
        from pynol.common.sequence.RNA import RNA

        all_features = []
        with open(file) as handle:
            full_gbk = [s for s in parse(handle)]
        for gbk in full_gbk:
            for feat in gbk.features:
                if not feat.key in {'source', 'gene'}:
                    if feat.key == "CDS":
                        all_features += [CDS.fromGBKGenomefile(feat, gbk)]
                    elif "RNA" in feat.key:
                        all_features += [RNA.fromGBKGenomefile(feat, gbk)]
                    else:
                        all_features += [Feature.fromGBKGenomefile(feat, gbk)]
        return all_features

    @classmethod
    def fromGBKGenomefile(cls, feature, gbk):
        obj = cls()
        loc = feature.location.replace('<','').replace('>','')
        rev_comp = "complement" in loc
        obj.strand = 1
        if rev_comp:
            loc = loc.replace("complement(",'').replace(')','')
            obj.strand = -1
        qualis =  { q.key[1:-1] : q.value.replace('"','') for q in feature.qualifiers}
        genomic = Genomic.find_one({'other_ids.original' : gbk.version})
        assert genomic, "the contigs for these features does not exists"
        start = int(loc.split('..')[0])
        end = int(loc.split('..')[1])
        seq = genomic.sequence[(start-1):end]
        extras = qualis
        extras['loc'] = feature.location
        seq.alphabet = DNAAlphabet()
        obj.type = feature.key
        obj.start = start
        obj.end = end
        obj.pretty_id = genomic.pretty_id+":" + obj.type + ":" + (extras['ID'].split("_")[-1] if extras.get('ID') else str(start) + "to" + str(end))
        obj.sequence = seq.reverse_complement() if rev_comp else seq
        obj.genomic = genomic.id
        obj.genome = genomic.genome
        obj.other_ids={}
        obj.more = extras
        obj.save()
        return obj

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
        obj.pretty_id = genomic.pretty_id+":" + obj.type + ":" + (extras['ID'].split("_")[-1] if extras.get('ID') else str(start) + "to" + str(end))
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
