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
from pandas import DataFrame, Index
from numpy import mean

class Blast( Tool, Slurm ):

    name = "Blast"
    executable = "blastn"
    version_cmd = executable + " -v"

    def run(self, sequence, reference, out_dir = "/tmp/", seq_typ = "prot") :
        exe_str = "{exe}  -num_threads {threads}  -query {sequnce}  -db {reference} -outfmt 10  >  {output}"
        threads = THREADS
        exe = "blastp" if  seq_typ == "prot" else "blastn"
        infile =  pjoin(out_dir, sequence.name + (".faa" if  seq_typ == "prot" else ".fna"))
        outfile = pjoin(out_dir, sequence.name + "_vs_" + os.path.basename(reference) + "." + exe)
        sequence.write_fasta(infile, pretty = True)
        return {'command' : exe_str.format(output = outfile, reference = reference, sequnce = infile, threads = threads, exe = exe), 'files' : [infile, reference]}

    def run_remote(self, sequence, reference, seq_typ = "prot"):
        stuff = self.run(sequence, sequence, '/home/moritz/temp/' + sequence.name, seq_typ)
        stuff['command'] = """module load bioinfo-tools
module load blast
mkdir {dones}
mkdir {path}
{command}
mv {path} {dones}
""".format(command = stuff['command'], path = '/home/moritz/temp/' + sequence.name, dones = '/home/moritz/temp/blasts/')

        Slurm.run_remote(self, job_name = "Blast_" + sequence.name, command = stuff['command'], files = stuff['files'], threads = 4)

    def retrieve_data(self, sequence, reference, seq_typ = "prot" ):
        exe = "blastp" if  seq_typ == "prot" else "blastn"

        file =  sequence.name + "_vs_" + os.path.basename(reference) + "." + exe
        Slurm.retrieve_data(self, files = [file], folder = '/home/moritz/temp/blasts/' + sequence.name)
        if not os.path.exists(pjoin("/tmp",file)):
            return "Not run or crashed"

    @classmethod
    def parse_blast_file(self, file, full = True, evalue_cutoff = 0.001, length_cutoff = 0.8):
        head =['query', 'subject', 'identity', 'length', 'mismatch', 'gaps', 'qstart', 'qend', 'sstart','send','evalue','bitscore']
        if os.stat(file).st_size > 0:
            raw_data = DataFrame.from_csv(file, header = -1, index_col = False)
        else :
            raw_data = DataFrame(columns=Index(head))
        raw_data.columns = Index(head)
        #m_length = mean([len(s) for s in self.get_seqs()])
        out_data = raw_data #.loc[raw_data['length'] > (m_length*length_cutoff)]
        out_data = out_data.loc[out_data['evalue'] < evalue_cutoff]
        return out_data
