# File: Sequence.py
from pynol.common.Record import Record
import hashlib
from Bio.Seq import Seq

class Sequence( Record ):

# Attributes: Instance
    def __init__(self):
        self.pretty_id = None
        self.other_ids = {}

    def checksum(self):
        return hashlib.md5(str(self.sequence).encode('utf-8')).hexdigest()

# Operations
