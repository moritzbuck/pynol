import os
from os.path import join as pjoin
import BTrees.OOBTree
from Bio.Alphabet import DNAAlphabet
from Bio.Seq import Seq
import ZODB, ZODB.FileStorage
import transaction
from pynol.common.sequence.Genomic import Genomic
from pynol.common.genome.Genome import Genome
from pynol.common.genome.sources.FromFile import FromFile
from pynol.common.genome.sources.FromRefSeq import FromRefSeq


base = pjoin(os.environ['HOME'], 'repos/moritz/pynol')

storage =        ZODB.FileStorage.FileStorage(pjoin(base,'test_data','mydata.fs'))
db = ZODB.DB(storage)
connection = db.open()
root = connection.root

if not root.__dict__['_root'].get('seqs') :
    root.seqs = BTrees.OOBTree.BTree()

from pynol.common.genome.Genome import Genome
from pynol.common.genome.sources.FromFile import FromFile
from pynol.common.genome.sources.FromRefSeq import FromRefSeq

test_genome = Genome()
test_genome.source = FromRefSeq("test_genome", 'GCA_000989175.1')
test_genome.populate()

transaction.commit()
