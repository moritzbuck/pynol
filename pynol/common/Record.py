import persistent
import uuid

class Record( persistent.Persistent ):

    """docstring for Record."""

    def __init__(self):
        super(Record, self).__init__()
        self.id = uuid.uuid4()
