import ZODB, ZODB.FileStorage
import transaction
from BTrees.OOBTree import BTree

class PynolDB(object):

    def __init__(self, storage_path):
        self.path = storage
        self.storage = ZODB.FileStorage.FileStorage(self.path)
        self.db = ZODB.DB(storage)
        self.connection = self.db.open()
        self.data_root = self.connection.root

        if not root.__dict__['_root'].get('Genome') :
            root.Genome = BTree()
        if not root.__dict__['_root'].get('Sequence') :
            root.Sequence = BTree()
        if not root.__dict__['_root'].get('Source') :
            root.Source = BTree()
