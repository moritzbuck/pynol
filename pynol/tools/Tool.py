import os
from mongo_thingy import Thingy
import distutils.spawn

class Tool( Thingy ):

    @classmethod
    def make(cls, executable = None):
        obj = cls()
        obj.name = cls.name
        obj.executable = executable if executable else cls.executable
        assert distutils.spawn.find_executable(obj.executable), "%s's executable %s not found" % (obj.name, obj.executable)
        return obj

    def run(self):
        pass
