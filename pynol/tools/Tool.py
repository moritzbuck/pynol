
import os
from mongo_thingy import Thingy
import distutils.spawn
import subprocess

class Tool( Thingy ):

    @classmethod
    def make(cls, executable = None):
        obj = cls()
        obj.name = cls.name
        obj.executable = executable if executable else cls.executable
        obj.version_cmd = cls.version_cmd
        assert distutils.spawn.find_executable(obj.executable), "%s's executable %s not found" % (obj.name, obj.executable)
        obj.local_version = obj.run_version_command()

        return obj

    def run(self):
        pass

    def run_version_command(self):
        return subprocess.getoutput(self.version_cmd)
