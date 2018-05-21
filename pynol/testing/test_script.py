    import os
    from os.path import join as pjoin
    from Bio.Alphabet import DNAAlphabet
    from Bio.Seq import Seq
    from pynol.common.sequence.Genomic import Genomic
    from pynol.common.genome.Genome import Genome
    from pynol.common.genome.sources.FromRefSeq import FromRefSeq
    from pynol.common.genome.sources.FromFile import FromFile

    from pynol.database.PynolDB import PynolDB
    import pynol
    from tqdm import tqdm
    from pynol.tools.gene_prediction.Prokka import Prokka


    from mongo_thingy import connect, Thingy
    connect("mongodb://localhost/pateci")

    with open("test_data/caulos.tsv") as handle:
        ids = [s.strip().replace('"', '').split("\t") for s in handle.readlines()]

    for s in tqdm(ids[1:]):
        if not "UBA" in s[0]:
            if not Genome.find_one({'other_ids.refseq_id' : s[0]}):
                print("Downloading ", s[0])
                genome = Genome()
                genome.name = "Ppi_" + s[0]
                genome.source = FromRefSeq.FromRefSeq(s[0])
                genome.taxonomy = {}
                genome.taxonomy['ncbi'] = s[2]
                genome.taxonomy['gtdb'] = s[3]
                genome.populate()
                genome.save()
        else :
            if not Genome.find_one({'other_ids.uba_id' : s[0]}):
                print("Catching ", s[0])
                genome = Genome()
                genome.name = "Ppi_" + s[0]
                pat = "/home/moritzbuck/uppmax/temp/release80/UBAs/{uba_id}.fsa".format(uba_id = s[0])
                genome.source = FromFile.FromFile(pat)
                genome.taxonomy = {}
                genome.taxonomy['ncbi'] = s[2]
                genome.taxonomy['gtdb'] = s[3]
                genome.other_ids = {}
                genome.other_ids['uba_id'] = s[0]
                genome.populate()
                genome.save()


def old():
    with open("test_data/patsci.ids") as handle:
         ids = [s.strip() for s in handle.readlines()]

    #remove because duplicate
    ids.remove('GCA_001776255.1')

    for s in tqdm(ids):
        if not Genome.find_one({'other_ids.refseq_id' : s}):
            genome = Genome()
            genome.name = "Patescibacteria_" + s
            genome.source = FromRefSeq.FromRefSeq(s)
            genome.populate()
            genome.save()

    with open("test_data/chitinoi.ids") as handle:
         ids = [s.strip() for s in handle.readlines()]

    for s in tqdm(ids):
        if not Genome.find_one({'other_ids.refseq_id' : s}):
            genome = Genome()
            genome.name = "Chitinos_" + s
            genome.source = FromRefSeq.FromRefSeq(s)
            genome.populate()
            genome.save()

    with open("test_data/comamos.ids") as handle:
        ids = [s.strip() for s in handle.readlines()]



    for s in tqdm(ids):
        if not Genome.find_one({'other_ids.refseq_id' : s}):
            genome = Genome()
            genome.name = "Comamo_" + s
            genome.source = FromRefSeq.FromRefSeq(s)
            genome.populate()
            genome.save()




tt = Genome.find_one()
pp = Prokka.make()
pp.run(tt)
