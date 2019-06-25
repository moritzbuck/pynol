# File: Genome.py
from pynol.common.sequence.Sequence import Sequence
from mongo_thingy import Thingy
import re
from bson.objectid import ObjectId
from pynol.common.sequence.sets.sources.FromFile import FromFile
import sys
from bson.objectid import ObjectId
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import hashlib
import os


class Alignment( Thingy ):


    def __len__(self):
        return len(self.seqs)

    def __getitem__(self, key):
        if type(key) == int :
            return self.seqs[key]

    def populate(self):
        self.source.get_data(self)


    def write_fasta(self, file, pretty = False):

        if pretty:
            to_seq_record = lambda c : SeqRecord(id = c.pretty_id, description = "", seq = c.seq)
        else :
            to_seq_record = lambda c : SeqRecord(id = str(c.id), description = "", seq = c.sequence)

        recs = [to_seq_record(c) for c in self.seqs]
        with open(file, "w") as handle:
            SeqIO.write(recs, handle, "fasta")

    def consensus(self):
        cons = ""
        for i in range(len(self.seqs[0])):
            aas = [s.sequence[i] for s in self.seqs]
            abv = {aa : aas.count(aa) for aa in set(aas)}
            cons += max(abv.items(), key = lambda a: a[1])[0]
        return cons

    def single_genome(self):
        if all([ObjectId(s.id) != None for s in  self.seqs]):
            genomes = {s.id for s in self.seqs}
            byGenomes = {g : [s for s in self.seqs if s.id == g] for g in genomes}
        else :
            genomes = {s.genome for s in self.seqs}
            byGenomes = {g : [s for s in self.seqs if s.genome == g] for g in genomes}
        pick_rng = lambda v : int(hashlib.md5(str(v[0].sequence).encode('utf-8')).hexdigest(), 16) % len(v)
        out_seqs = [SeqRecord(id = str(k), seq=v[pick_rng(v)].sequence, description="", name = "") for k, v in byGenomes.items()]
        return out_seqs



    @property
    def seqs(self):
        if not self._seqlist:
            self._seqlist = list(Sequence.find({'seqset' : self.id}))
        return self._seqlist

    def block_seqs(self, fract_cols = 0.5, fract_rows = 0.5):
        mask = [ (1-[s.sequence[i]   for s in self.seqs].count('-')/len(self.seqs)) > 0.5 for i in range(len(self.seqs[0])) ]
        sts = [ SeqRecord(id = str(s.ori_cds) if s.ori_cds else str(s.id), description = "", seq = Seq("".join([ s.sequence[i] for i, m in zip(range(len(self.seqs[0])), mask) if m ])) )for s in self.seqs]
        sts = [s for s in sts if (1-str(s.seq).count("-")/len(s.seq)) > fract_rows]
        return sts

    @classmethod
    def cat(cls, alignments, genome_ids = None):
        alis = [{ str(Sequence.find_one(ObjectId(s.id)).genome) : s for s in ali.single_genome() } for ali in alignments]
        ali_lens = [len(ali.seqs[0]) for ali in alignments]

        if not genome_ids :
            genome_ids = set(sum([ list(ali.keys()) for ali in alis], []))

        out_seqs = { }
        for g_id in genome_ids:
            s = [str(ali.get(g_id).seq) if ali.get(g_id) else '-'*ali_len for ali, ali_len in zip(alis, ali_lens)]
            out_seqs[g_id] = "".join(s)

        return [SeqRecord(id = str(k) , seq = Seq(s), description = "", name = "")for k,s in out_seqs.items()]


    @classmethod
    def FromRawFile(cls,name, file_name, whatever_id, additional_info = {}):
        alignmnt = cls()
        alignmnt.save()

        alignmnt.name = name
        alignmnt.source = FromFile.FromFile(file_name)
        alignmnt.save()
        print("Created alignmnt with {id}".format(id = alignmnt.id))
        alignmnt.populate()
        return alignmnt


    @classmethod
    def FromSeqRecs(cls,parent, seqs):
        alignmnt = cls()
        alignmnt.save()
        alignmnt.parent = [p.id for p in parent]
        alignmnt.name = ":".join([p.name for p in parent])
        SeqIO.write(seqs, "1234xcd.fasta", "fasta")
        alignmnt.source = FromFile.FromFile("1234xcd.fasta")
        alignmnt.save()
        print("Created alignmnt with {id}".format(id = alignmnt.id))
        alignmnt.populate()
        os.remove("1234xcd.fasta")
        return alignmnt




    @property
    def source(self):
        source_class = getattr(sys.modules[__name__], self._source_class)
        return source_class.find_one({ '_id' : self._source})

    @source.setter
    def source(self, source):
        source.save()
        self._source = source._id
        self._source_class = source.__class__.__name__
