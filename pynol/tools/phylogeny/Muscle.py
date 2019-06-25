from pynol.tools.Tool import Tool
from pynol.common.genome.Genome import Genome
import uuid
from pynol import THREADS
from pynol.tools.slurm import Slurm
from os.path import join as pjoin
from Bio import SeqIO
import os
from pynol.common.sequence.Feature import Feature
from pynol.common.sequence.RNA import RNA
from pynol.common.sequence.CDS import CDS
from pynol.common.sequence.Genomic import Genomic
from tqdm import tqdm
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

class Muscle( Tool ):

    name = "Muscle"
    executable = "muscle"
    version_cmd = executable + " -version"

    def run(self, sequences, seq_names,  out_dir = "/tmp/", name = None) :

        exe_str = "muscle -out {outfile}  -in {infile} "

        if not name:
            name = str(uuid.uuid4()).split("-")[-1]

        infile =  pjoin(out_dir, name +  ".fasta")
        outfile =  infile.replace(".fasta", "_aligned.fasta")
        SeqIO.write([SeqRecord(id = sname, seq = Seq(str(seq)), description="") for seq, sname in zip(sequences, seq_names)], infile, "fasta")

        os.system(exe_str.format(outfile = outfile, infile = infile))

        return [s for s in SeqIO.parse(outfile, "fasta")]
