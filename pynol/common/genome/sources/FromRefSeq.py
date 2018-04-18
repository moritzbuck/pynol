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

    def __init__(self, name, refseq_id, long_name = None):
        super(FromRefSeq, self).__init__(name, long_name)
        self.refseq_id = refseq_id


    def download_assembly(self):
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
        assert len(ass_id) == 1
        ass_id = ass_id[0]
        ftp.cwd(ass_id)
        all_files = ftp.nlst()
        assembly_file = [f for f in all_files if f.endswith(".fna.gz") and not "_cds_" in f]
        assert len(assembly_file) == 1
        assembly_file = assembly_file[0]
        stream = io.BytesIO()
        ftp.retrbinary('RETR %s' % assembly_file, stream.write)
        stream.seek(0)
        stream = FakeFile(gzip.GzipFile(fileobj=stream, mode = "rb"))
        return SeqIO.parse(stream, "fasta", alphabet = DNAAlphabet())

    def get_data(self):
        contigs = []
        for s in self.download_assembly():
            c = Genomic()
            c.sequence = s.seq
            c.type = "contig"
            c.other_ids['original'] = s.id
            contigs += [c]

        pads = len(str(len(contigs)))
        for i, s in enumerate(contigs):
            c.pretty_id = "{name}:contig:{i}/{len}".format(name = self.name, i = str(i+1).zfill(pads), len = len(contigs))

        return contigs
