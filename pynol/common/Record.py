import persistent
import uuid
import pynol

class Record( persistent.Persistent ):

    """docstring for Record."""

    def __init__(self, database = None):
        super(Record, self).__init__()
        self.id = str(uuid.uuid4())
        if database :
            self.db = database
        else :
            self.db = pynol.PYNOL
        self.db.full[self.id] = self
