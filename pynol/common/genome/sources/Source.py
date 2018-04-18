# File: Source.py
from pynol.common.Record import Record
from abc import ABC, abstractmethod
import re

class Source( Record, ABC ):

# Attributes: Instance

    def __init__(self, name, long_name):
        super( Record , self).__init__()
        self.name = name
        self.long_name = long_name
# Operations

    @abstractmethod
    def get_data(self):
        pass

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        assert len(value) < 36, "name is too long, use 'long_name' for this"
        assert bool(re.match('^[A-Za-z0-9_.]+$', value)), "please use only alphanumeric characters of underscore"
        self._name = value
