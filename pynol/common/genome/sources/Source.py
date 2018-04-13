# File: Source.py
from pynol.common.Record import Record
from abc import ABC, abstractmethod

class Source( Record, ABC ):

# Attributes: Instance

    def __init__(self):
        super( Record , self).__init__()
# Operations

    @abstractmethod
    def get_data(self):
        pass
