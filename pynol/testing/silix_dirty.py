import uuid
from os.path import join as pjoin

wk_path = pjoin("/tmp/",str(uuid.uuid4()))
os.makedirs(wk_path)

genomes = tt
for g in tqdm(genomes) :
    g.write_fasta(pjoin(wk_path, str(g.id) + ".faa"))


fasta_path = pjoin(wk_path, "all_proteoms.faa")
bast_path = pjoin(wk_path, "all_proteoms.more_sensitive.dblastp")
fastas = [pjoin(wk_path,f) for f in os.listdir(wk_path) if ".faa" in f and not "full" in f]
with open(fasta_path, 'w') as outfile:
    for fname in tqdm(fastas):
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)


os.system("diamond makedb --in {file} --db {file}".format(file = fasta_path))

os.system("diamond  blastp  --more-sensitive -p 20 -f  6 -q {file} --db {file} -o {out}".format(file =  , out = blast_path))

os.system("silix {fasta} {blast}".format(fasta = fasta_path, blast = blast_path))
