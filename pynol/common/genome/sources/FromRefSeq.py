from pynol.common.genome.sources.Source import Source
from pynol.common.sequence.Genomic import Genomic
from Bio import SeqIO
from Bio.Alphabet import DNAAlphabet
import re
import sys
import ftplib
import os
from ftplib import FTP
import io
import gzip
from io import StringIO

class FakeFile():
    def __init__(self, gzip_stream):
        self.stream = gzip_stream

    def readline(self):
        return self.stream.readline().decode()


class FromRefSeq(Source):
    """docstring for FromRefSeq."""

    @classmethod
    def FromRefSeq(cls,refseq_id):
        obj = cls()
        obj.refseq_id = refseq_id
        obj.save()
        return obj

    @property
    def refseq_id(self):
        return self._refseq_id

    @refseq_id.setter
    def refseq_id(self, ref_id):
        self._refseq_id = ref_id

    def get_seqio_genome_parser(self, type = "assembly"):
        ftp_addr = "ftp.ncbi.nlm.nih.gov"
        ftp_path = "genomes/all/{letters}/{first}/{second}/{third}/"

        letters = self.refseq_id.split("_")[0]
        numbers =  self.refseq_id.split("_")[1].split(".")[0]
        first, second, third = numbers[0:3], numbers[3:6], numbers[6:9]
        ftp_path = ftp_path.format(letters = letters, first=first, second=second, third=third)

        ftp=FTP(ftp_addr)
        ftp.login()
        ftp.cwd(ftp_path)
        ass_id = ftp.nlst()
        ass_id = [a for a in ass_id if self.refseq_id in a]
        assert len(ass_id) == 1
        ass_id = ass_id[0]
        ftp.cwd(ass_id)
        all_files = ftp.nlst()
        if type == "assembly" :
            assembly_file = [f for f in all_files if f.endswith(".fna.gz") and not "_cds_" in f and not "_rna_" in f]
        elif type == "gff" :
            assembly_file = [f for f in all_files if f.endswith("_genomic.gff.gz")]
        else :
            assert False, "Bad type"

        assert len(assembly_file) == 1
        assembly_file = assembly_file[0]
        stream = io.BytesIO()
        ftp.retrbinary('RETR %s' % assembly_file, stream.write)
        stream.seek(0)
        stream = FakeFile(gzip.GzipFile(fileobj=stream, mode = "rb"))
        if type == "assembly" :
            return SeqIO.parse(stream, "fasta", alphabet = DNAAlphabet())
        elif type == "gff" :
            return SeqIO.parse(stream, "gff", alphabet = DNAAlphabet())
        if not genome.other_ids:
            genome.other_ids = {}
        genome.other_ids['refseq_id'] = self.refseq_id
