    import os
    from os.path import join as pjoin
    from Bio.Alphabet import DNAAlphabet
    from Bio.Seq import Seq
    from pynol.common.sequence.Genomic import Genomic
    from pynol.common.genome.Genome import Genome
    from pynol.common.genome.sources.FromRefSeq import FromRefSeq
    from pynol.common.genome.sources.FromFile import FromFile
    from pynol.common.taxonomy.Taxonomy import Taxonomy
    from pynol.database.PynolDB import PynolDB
    import pynol
    from tqdm import tqdm
    from pynol.tools.gene_prediction.Prokka import Prokka


    from mongo_thingy import connect, Thingy
    connect("mongodb://localhost/pateci")

    taxonomy = Taxonomy()
    taxonomy.from_gtdb_taxonomy_file("/home/moritzbuck/data/bac_taxonomy_r80.tsv")

    test_genome = Genome.find_one({'_name' : "UBA9650"})

    prok = Prokka.make()

    with open("test_data/caulos.tsv") as handle:
        ids = [s.strip().replace('"', '').split("\t") for s in handle.readlines()]

    for s in tqdm(ids[1:]):
        if not "UBA" in s[0]:
            if not Genome.find_one({'other_ids.refseq_id' : s[0]}):
                print("Downloading ", s[0])
                taxon_str = s[3]
                assert taxonomy.check_consistency(taxon_str)
                Genome.FromRefSeq(s[0], taxon_str, s[2])

        else :
            if not Genome.find_one({'other_ids.uba_id' : s[0]}):
                print("Catching ", s[0])
                taxon_str = s[3]
                taxonomy.check_consistency(taxon_str)
                pat = "/home/moritzbuck/uppmax/temp/release80/UBAs/{uba_id}.fsa".format(uba_id = s[0])
                genome = Genome.FromUBAFile(pat, s[0], taxon_str)

    search_string = "^{tax_string}.*".format(tax_string = taxonomy['p__Patescibacteria'].get_tax_string(full = True))
    Genome.find({'taxonomy.gtdb' : {'$regex' : search_string}})

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
