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
mkdir {dones}
mkdir {path}
{command}
mv {path} {dones}
""".format(command = stuff['command'], path = '/home/moritz/temp/' + genome.name, dones = '/home/moritz/temp/prokkas/')

        Slurm.run_remote(self, job_name = "Prokka_" + genome.name, command = stuff['command'], files = stuff['files'], threads = 4)

    def retrieve_data(self, genome):
        gbk_file = pjoin('/tmp/', genome.name + ".gff")
        Slurm.retrieve_data(self, files = [genome.name + ".gff"], folder = '/home/moritz/temp/prokkas/' + genome.name)
        if not os.path.exists(gbk_file):
            return "Not run or crashed"
        with open(gbk_file) as handle:
            lines = [l.strip().split("\t") for l in handle if l[0] != "#" and "\t" in l
            ]

        all_feats = {'CDS' : {}, 'RNA' : {}, 'other' : {}}

        for l in lines :
            if l[2] == "CDS" :
                feat = CDS.fromGFFline(l)
                all_feats['CDS'][feat.pretty_id] = feat.id
            elif "RNA" in l[2] :
                feat = RNA.fromGFFline(l)
                all_feats['RNA'][feat.pretty_id] = feat.id
            else :
                print("Feature of type ", l[2], " with no specific class")
                feat = Feature.fromGFFline(l)
                all_feats['other'][feat.pretty_id] = feat.id

        return all_feats
