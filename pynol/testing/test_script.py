import os
from os.path import join as pjoin
import BTrees.OOBTree
from Bio.Alphabet import DNAAlphabet
from Bio.Seq import Seq
import ZODB, ZODB.FileStorage
import transaction
from pynol.common.sequence.Genomic import Genomic
from pynol.common.genome.Genome import Genome
from pynol.common.genome.source import FromFile


base = pjoin(os.environ['HOME'], 'repos/moritz/pynol')

storage =        ZODB.FileStorage.FileStorage(pjoin(base,'test_data','mydata.fs'))
db = ZODB.DB(storage)
connection = db.open()
root = connection.root

if not root.__dict__['_root'].get('seqs') :
    root.seqs = BTrees.OOBTree.BTree()

bla = Genomic()
bla.sequence = Seq("ATGCGCGCG", DNAAlphabet())
transaction.commit()
