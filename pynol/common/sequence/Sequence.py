# File: Sequence.py
import persistent

class Sequence( persistent.Persistent ):
    count = 0

# Attributes: Instance


    def __init__(self):
        self.id = Sequence.count
        Sequence.count += 1
        self.pretty_id = None
        self.other_ids = {}

# Operations
