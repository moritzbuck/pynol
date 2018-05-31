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


class Prokka( Tool, Slurm ):

    name = "SILIX"
    executable = "silix"
    version_cmd = executable + " -v"

    def run(self, genome, out_dir = "/tmp/") :
        exe_str = "prokka --outdir {out_dir}  --force --prefix {prefix} --locustag {tag} --cpus {threads} {infile} "
        prefix = genome.name
        tag = prefix
        threads = THREADS
        infile =  genome.name + ".fasta"
        genome.to_fasta(infile)
        return {'command' : exe_str.format(out_dir = out_dir, prefix = prefix, tag = tag, threads = threads, infile = infile), 'files' : [infile]}
