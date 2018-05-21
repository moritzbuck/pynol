from pynol.tools.Tool import Tool
from pynol.common.genome.Genome import Genome
import uuid
from pynol import THREADS
from pynol.tools.slurm import Slurm
from os.path import join as pjoin
from Bio import SeqIO
import os

class Prokka( Tool, Slurm ):

    name = "Prokka"
    executable = "prokka"

    def run(self, genome, out_dir = "/tmp/") :
        exe_str = "prokka --outdir {out_dir}  --force --prefix {prefix} --locustag {tag} --cpus {threads} {infile} "
        prefix = genome.name
        tag = prefix
        threads = THREADS
        infile =  genome.name + ".fasta"
        genome.to_fasta(infile)
        return {'command' : exe_str.format(out_dir = out_dir, prefix = prefix, tag = tag, threads = threads, infile = infile), 'files' : [infile]}

    def run_remote(self, genome):
        stuff = self.run(genome, '/home/moritz/temp/' + genome.name)
        stuff['command'] = """module load bioinfo-tools
module load prokka
mkdir {path}
{command}""".format(command = stuff['command'], path = '/home/moritz/temp/' + genome.name)

        Slurm.run_remote(self, job_name = "Prokka_" + genome.name, command = stuff['command'], files = stuff['files'], threads = 4)

    def retrieve_data(self, genome):
        gbk_file = pjoin('/tmp/',genome.name, genome.name + ".gbk")
        Slurm.retrieve_data(self, files = [genome.name], folder = '/home/moritz/temp/')

        #clean prokka but with gbk header
        nb_ctgs = len(genome.contigs)

        sed = "sed -i 's/\/{counts}/\/{counts}\t/' {gbk}"

        sed = sed.format(counts = nb_ctgs, gbk = gbk_file)

        os.system(sed)

        return [s for s in SeqIO.parse(gbk_file, "genbank")]
